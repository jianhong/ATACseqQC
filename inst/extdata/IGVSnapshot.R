#' Take IGV snapshots
#' @description Take IGV snapshots for genomic regions of interest.
#' Prerequisite: Make sure the Integrative Genomics Viewer (IGV) 
#' (http://software.broadinstitute.org/software/igv/) is installed in a 
#' user's computer. Be aware of the MAC OSX security and privacy setting. 
#' To utilize this functionality, be sure that the box for "Enable port", which 
#' can be found under the toolbar (View -> Preference -> Advanced), is checked 
#' to allow communication via port.
#' @param maxMem A character vector of length 1L containing the maximum usable 
#' memory for the IGV to be launched, which is defined as the following: 
#' 'mm' - 1.2 GB , 'lm' - 2 GB, 'hm' - 10 GB, ' - 750 MB. For details, please
#' see the #' parameter memory of function startIGV in the SRAdb package.
#' @param genomeBuild A character vector of length 1L: the genome build version.
#' For example, "hg38" for the human reference genome assembly GRCh38. 
#' Please see IGV manual for supported genomes of species. 
#' @param bamFileFullPathOrURLs A character vector containing the full paths 
#' or ULRs for bam files and matching bam index files. No spaces are allowed 
#' in the paths or URLs.
#' @param geneNames A character vector containing the names of genes or genomic
#' regions for visualization. For human (maybe other mammalian) ATAC-seq data, 
#' names of several house-keeping genes are provided as defaults.
#' @param genomicRegions A character vector containing genomic regions in the
#' form of "chr1:20000-40000".
#' @param sessionFile A character vector of length 1L containing the name of 
#' the IGV session file to be saved.
#' @param outDir A character vector of length 1L containing the output directory.
#' @author Haibo Liu
#' @return None
#' @export
#' @importFrom utils installed.packages
#' @details This function depend on SRAdb startIGV IGVsocket IGVgenome IGVload IGVgoto IGVsnapshot IGVsession
#' @examples
#' if(interactive()){
#' library(SRAdb)
#' exampleBams = file.path(system.file('extdata',package='SRAdb'),
#'              dir(system.file('extdata',package='SRAdb'),pattern='bam$'))
#' IGVSnapshot(maxMem="mm", genomeBuild="hg18", bamFileFullPathOrURLs= exampleBams,
#'      geneNames=NULL, genomicRegions="chr1:1-1000",  sessionFile="IGV.session")
#'}


IGVSnapshot <- function(maxMem, genomeBuild="hg38", bamFileFullPathOrURLs, 
                        geneNames=c("GAPDH", "ACTB", "C1orf43", "CHMP2A", "EMC7", 
                                    "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", 
                                    "VCP", "VPS29"), genomicRegions, sessionFile, outDir=getwd())
{
  if(!"SRAdb" %in% installed.packages()[, 1]){
    stop("SRAdb is required for IGVSnapshot.")
  }
    SRAdb::startIGV(maxMem)
    while (TRUE && interactive())
    {
        setupIGV <- readline(prompt="Is IGV starting up and port enabled? Y or N \n")
        if (grepl("Y|y", setupIGV, ignore.case=TRUE, perl=TRUE ))
        {
            break
        }else{
            print("Please set up IGV and enable port! For MAC OSX, please check the
                  security and privacy setting to allow IGV to start up regardless 
                  of the security setting. Please see the IGV manual for details.
                  Please also adjust the IGV window to full screen.")
        }
    }
    
    sock <- SRAdb::IGVsocket()
    SRAdb::IGVgenome(sock, genome=genomeBuild)
    SRAdb::IGVload(sock, bamFileFullPathOrURLs)
    
    if (!missing(genomicRegions))
    {
        genomicRegions <- genomicRegions[grepl("^\\w+:\\d+-\\d+$", genomicRegions, perl=TRUE)]
        geneNames <- c(geneNames, genomicRegions)
    }
    
    if (length(geneNames)>=1)  ## there is at least one gene name or one genomic region in correct forms
    {
        ## Take snapshot for interesting genomic regions
        for (gene in geneNames)
        {
            takeSnapshot(sock=sock, gene=gene, outDir=outDir)
        }
    }else if (interactive()){  ## ask to the user to input gene names or genomic regions
        takeSnapshotOUponInput(sock, outDir)
    }else{
        stop("At least one geneName or genomicRegion should be provided!")
    }
    
    ## save the IGV session file
    if (interactive())
    {
        saveSession <- readline(prompt="Save the IGV sessionfile? Y or N \n")
        if (grepl("Y|y", saveSession, perl=TRUE))
        {
          SRAdb::IGVsession(files=bamFileFullPathOrURLs, sessionFile=sessionFile, 
                       genome=genomeBuild, VisibleAttribute='', destdir=outDir)
        }
    }else{
      SRAdb::IGVsession(files=bamFileFullPathOrURLs, sessionFile=sessionFile, 
                   genome=genomeBuild, VisibleAttribute='', destdir=outDir)
    }
    
}

#### Helper function: interactively get geneNames or genomicRegions from standard input, 
#### and take a snapshot after asking.
takeSnapshotOUponInput <- function(sock, gene, outDir)
{
    while (TRUE)
    {
        input <- readline(prompt="Please input gene a name (such as ACTB), genomic region 
                          (such as chr1:10000-20000) or type 'Y' to quit: \n")
        if (grepl("\\w+|\\w+:\\d+-\\d+", input, perl=TRUE))
        {
            takeSnapshot(sock=sock, gene=gene, outDir=outDir)
        }else if (grepl("Y|y", input, ignore.case=TRUE, perl=TRUE)){
            break
        }else{
            readline(prompt="The input is not valid. Please input gene a name (such as 'ACTB'), 
                     genomic region (such as 'chr1:10000-20000') or type 'Y' to quit: \n")
        }
        }
    }


#### helper function to take a IGV snapshot after asking
takeSnapshot <- function(sock, gene, outDir)
{
  SRAdb::IGVgoto(sock, gene)
    if (interactive())
    {
        takeSnapshot <- readline(prompt="Take a snapshot? Y or N \n")
        if (grepl("Y|y", takeSnapshot, perl =TRUE))
        {
          SRAdb::IGVsnapshot(sock)
        }
    }else{
      SRAdb::IGVsnapshot(sock)
    }
}

