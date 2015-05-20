#!/usr/bin/RScript
# htseq-combine.R
# copy this code to RStudio and adapt file locations to match yours

# Take 'gtf' htseq-count results and melt them in to one big dataframe

# where are we?
basedir <- "/work/Projects-BITS/PJ1502_WP_comparing-mappers/GIB-counts"
setwd(basedir)

cntdir <- paste(basedir, "htseq_counts", sep="/")
pat <- "-counts.txt"
htseq-counts <- list.files(path = cntdir,
	pattern = pat,
	gtf.files = TRUE,
	recursive = FALSE,
	ignore.case = FALSE,
	include.dirs = FALSE)

myfiles <- htseq-counts
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
	infile = paste(cntdir, myfiles[i], sep = "/")
	DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
	cnts <- gsub("(.*)_gtf_counts.txt", "\\1", myfiles[i])
	colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge gtf elements based on first ID columns
data <- DT[[myfiles[1]]]

# inspect
head(data)

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
	y <- DT[[myfiles[i]]]
	z <- merge(data, y, by = c("ID"))
	data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# remove no-count rows
data <- data[data$tot.counts>0,]

# inspect and look at the top row names!
head(data)
tail(data)

####################################
# take summary rows to a new table
# ( not starting with ENS with invert=TRUE )

# transpose table for readability
data.summary <- data[grep("^ENS", rownames(data), perl=TRUE, invert=TRUE), ]

# review
data.summary

# transpose table
t(data.summary)

# write summary to file
write.csv(data.summary, file = "htseq_counts-summary.csv")

####################################
# take gtf data rows to a new table

data.table <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# final merged table
head(data.table, 3)

# write data to file
write.csv(data.table, file = "htseq_counts.csv")

# cleanup intermediate objects
rm(y, z, i, DT)
