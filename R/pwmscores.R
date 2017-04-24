#' @title max PWM scores for sequences
#' @description calculate the maximal PWM scores for each given sequences
#' @param pwm A Position Weight Matrix represented as a numeric matrix
#'            with row names A, C, G and T.
#' @param subject Typically a \link[Biostrings]{DNAString} object.
#'            A \link[IRanges]{Views} objecton
#'            a \link[Biostrings]{DNAString} subject,
#'            a \link[Biostrings]{MaskedDNAString} object,
#'            or a single character string, are also supported.
#'            IUPAC ambiguity letters in \code{subject} are ignored (i.e.
#'            assigned weight 0) with a warning.
#' @importFrom Biostrings DNAString PWMscoreStartingAt
#' @return a numeric vector
#' @author Jianhong

pwmscores <- function(pwm, subject){
  len <- length(subject)
  subject <- split(subject,
                   rep(seq_len(ceiling(len/1e5)), each=1e5)[seq_len(len)])
  motif.len <- ncol(pwm)
  score <- lapply(subject, function(.subject){
    .len <- length(.subject)
    .subject <- DNAString(paste(as.character(.subject), collapse=""))
    s <- PWMscoreStartingAt(pwm = pwm, subject = .subject,
                            starting.at = seq.int(80 * .len - motif.len))
    s <- split(s, rep(seq_len(.len), each=80)[seq_along(s)])
    s <- lapply(s, function(.ele) max(.ele[seq_len(80-motif.len)]))
    unname(unlist(s[order(as.numeric(names(s)))]))
  })
  unname(unlist(score[order(as.numeric(names(score)))]))
}
