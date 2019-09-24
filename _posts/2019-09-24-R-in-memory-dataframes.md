# data.frames on steroids

### Some background

Some days ago landed on my desk a request to identify bidirectional promoters and plot the signal of an NGS experiment around these promoters.

We were interested in identifying gene pairs following one of these conformations (classes):

```
A) head-head (A1: divergent, A2: convergent, A3: same strand)
  1)         +--->        2)       +--->            3)   +--->
   ----------+----------   --------+-------------        | +--->
   --------+------------   ----------+-----------      --+-+------ ... -----------
       <---+                     <---+                 ----------- ... ------+-+--
                                                                         <---+ |
                                                                           <---+
B) tail-head
    +---> +--->
   -+-----+----- ... -------------
   ------------- ... -----+-----+-
                      <---+ <---+
C) tail-tail
      +--->
   ---+-----------
   -----------+---
          <---+
```

One can certainly find lists with bidirectional promoters as supplementary of some manuscripts, but I had some reservations in using those lists:
  * some were outdated, or not available anymore
  * I was interested in the most common isoform of a gene only, and wanted to use a modern annotation of that (see [APPRIS](http://appris.bioinfo.cnio.es/) and [MANE Select](https://ncbiinsights.ncbi.nlm.nih.gov/2019/03/12/mane-select-v0-5/))
  * I needed custom types of bidirectionality, involving also the 3' end (see classes B and C)

For these reasons I needed to calculate them from scratch.

### Getting the list of transcripts

Getting the most abundant isoform is worth another post (and out of scope here). Oversimplifying a bit, get from BiomaRt a table with gene and transcript name, and genomic coordinates of all genes. Then filter out those non-PRINCIPAL isoforms annotated in APPRIS.

### Getting the list of bidirectional pairs of genes
There's several ways one can code something like this. The question is how to do this quickly, in terms of coding and execution time.

First though was using the fantastic [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package from Bioconductor. It has plenty of methods to work with genomic coordinates. The main drawback is that I'm not very fluent with GenomicRanges, and after gathering all my knowledge and going throught the vignette, I didn't identify a simple way to do this. Same with the *tidyverse* verbs, as you may suggest.

I won't extend any further regarding other possible ways to do this (BEDtools?). Just on the aim of this post: R data.frames (plus a bit of help from SQL).

Keep in mind the starting point is tabulated data from BiomaRt, that's been crossed and filtered with data from APPRIS. This part is done within R, therefore we now have a data.frame containing all principal isoforms + coordinates. We can start writing some coordinate arithmetic within `[]` to get the bidirectional pairs, but that would be *insanely slow*. Also, we could try to figure out how to do it with *GenomicRanges*, the *tidyverse* or export the table from R and process it somewhere else (like *BEDtools* or *SQLite*).

Hmmm, *SQLite*! Writing such query in a relational database would be trivial, and the database backend fast as lightning. The idea of the cartesian product behind relational databases (each row in the first table is paired with all the rows in the second table) perfectly matches our problem. We can even have indexes which would speedup the search if needed.

Now, let's do it all without exiting R, nor exporting anything out of it. Simply *upgrading* our good old data.frame into a fully fledged table in an perishable (in-memory) database.

#### Geting stuff from BiomaRt and APPRIS

Let's get all principal isoforms from BiomaRt and APPRIS:

```R
appris <- read.delim("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt", head=FALSE)
regions  <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "chromosome_name",
                               "transcript_gencode_basic", "transcript_start", "transcript_end", "strand"),
                  filters="ensembl_transcript_id",
                  values=appris$V3[appris$V5 == "PRINCIPAL:1"],
                  mart=useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org")))
regions$tss <- ifelse(regions$strand > 0, regions$transcript_start, regions$transcript_end)  # Biomart encodes '+' as 1 and '-' as -1
```

#### Write this stuff to an in-memory db, to do quick and simple cartesian products

```R
library(DBI)
db <- dbConnect(RSQLite::SQLite(), ":memory:")  # create an ephemeral in-memory RSQLite database
dbWriteTable(db, "regions", regions)            # create a table with the same structure and contents as in the data.frame
```

#### Get the bidirectional pairs
The SQL query is damned simple and as explicit as English: get all gene pairs on the same chromosome, different gene\_id, different strand, which the distance between TSSes is less than 1000bp.

```R
query <- paste("SELECT *",
               "  FROM regions AS gene1, regions AS gene2",
               " WHERE gene1.chromosome_name = gene2.chromosome_name",
               "   AND gene1.ensembl_gene_id <> gene2.ensembl_gene_id",
               "   AND gene1.strand <> gene2.strand",
               "   AND ABS(gene1.tss - gene2.tss) < 1000")

bidirectional <- dbGetQuery(db, query)
colnames(bidirectional) <- c(paste0(colnames(regions), ".gene1"),
                             paste0(colnames(regions), ".gene2"))
write.csv(bidirectional, file="bidirectional_promoters.csv", row.names=FALSE)
```

*Et voila!* less than 10 seconds to run.

### Conclusions

*in-memory* relational databases are great for doing some complex operations in basic data.frames (tabulated with simple datatypes), which otherwise would be slow or impossible to write. They're also simple and quick to construct within R from a data.frame. SQL is also a standard query language, turing complete and well understood. Therefore, data.frames converted into relational tables represent a powerful tool for many common problems in genomics.

### Addendum

No indexes are created with `dbWriteTable()`. That could have a noticeable performance issues with very big data.frames or intensive queries.

For instance:

```R
library(DBI)
db <- dbConnect(RSQLite::SQLite(), ":memory:")
dbWriteTable(db, "mtcars", mtcars)
head(dbReadTable(db, "mtcars"))
#    mpg cyl disp  hp drat    wt  qsec vs am gear carb
# 1 21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
# 2 21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
# 3 22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
# 4 21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
# 5 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
# 6 18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
```

To get the *path* that the DBMS will follow with a simple `WHERE` clause:

```R
dbGetQuery(db, "EXPLAIN QUERY PLAN SELECT COUNT(*) FROM mtcars WHERE mpg < 20")
#   id parent notused            detail
# 1  3      0       0 SCAN TABLE mtcars   <-- it's going to crawl through whole the table!
```

*SCAN TABLE is bad*.

So let's create an index on the `WHERE` group:

```R
dbExecute(db, "CREATE INDEX idx_mpg ON mtcars(mpg)")
# [1] 0
dbGetQuery(db, "EXPLAIN QUERY PLAN SELECT COUNT(*) FROM mtcars WHERE mpg < 20")
#   id parent notused                                                    detail
# 1  3      0       0 SEARCH TABLE mtcars USING COVERING INDEX idx_mpg (mpg<?)
```

*SEARCH TABLE is good*. This means it avoids going through the whole table, just the rows it knows mpg is below 20 (extrapolate it to a table with milions of rows where only few rows have mpg below 20).
