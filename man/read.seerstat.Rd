\name{read.seerstat}
\alias{read.seerstat}
\title{
Read data from SEER*Stat export files
}
\description{
A function to read data and dictionary information from SEER*Stat export files. 
It should be called to read in group survival data 
before running function \code{recurrisk.group}.
}
\usage{
read.seerstat(file.dic, file.txt)
}
\arguments{
  \item{file.dic}{
   The dictionary file exported from SEER*Stat software with .dic extension 
   which contains the information describing the layout of the export data file.
}
  \item{file.txt}{
   The cause-specific or relative survival data generated from SEER*Stat in .txt format.Note that, if the .txt file is not exported from SEER*Stat
   but manually converted, the user should make sure the delimiter used in the txt file is the same as the one specified in the dic file - field delimiter option. 
}
}
\details{
This function makes use of function read.SeerStat in R package SEER2R. 
}
\value{
A data frame containing SEER*Stat survival data.
}

\seealso{
\code{\link{recurrisk.group}}
}
%%\examples{
%%data.group<-read.seerstat("groupdata_example_seer.dic","groupdata_example_seer.txt")
%%}
