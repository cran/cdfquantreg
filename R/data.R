#' Ambiguity-Conflict data
#'
#' A data from a study that investigates the judgment under ambiguity and conflict 
#'
#' @format A data frame with 166 rows and 2 variables:
#' \describe{
#'   \item{ID}{subject ID}
#'   \item{value}{Rating in each judgment scenario}
#'   \item{scenario}{Index for judgment scenarios}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/16594767/}
"Ambdata"

#' Stress-Anxiety data
#'
#' A data from a study that investigates the relationship between stress and anxiety.
#'
#' @format A data frame with 166 rows and 2 variables:
#' \describe{
#'   \item{Anxiety}{Scores on Anxiety subscale}
#'   \item{Stress}{Scores on Stress subscale}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/16594767/}
"AnxStrData"


#' Juror data
#'
#' Juror Judgment Study.
#'
#' @format A data frame with 104 rows and 3 variables:
#' \describe{
#'   \item{crc99}{The ratings of confidence levels with rescaling into the (0, 1) interval to avoide 1 and 0 values.}
#'   \item{vert}{ was the dummy variable for coding the conditions of verdict types, whereas }
#'   \item{confl}{ was the dummy variable for coding the conflict conditions}
#' }
#' @source \doi{10.1375/pplt.2004.11.1.154}
"JurorData"


#' IPCC data-set
#'
#' The IPCC data-set comprises the lower, best, and upper estimates 
#' for the phrases "likely" and "unlikely" in six IPCC report sentences.
#'
#' @format A data frame with 4014 rows and 8 variables:
#' \describe{
#'   \item{subj}{Subject ID number}
#'   \item{treat}{Experimental conditions}
#'   \item{valence}{Valence of the sentences}
#'   \item{prob}{raw probability estimates} 
#'   \item{probm}{Linear transformed prob into (0, 1) interval} 
#'   \item{mid}{Distinguish lower, best and upper estiamtes }
#'   \item{high}{Distinguish lower, best and upper estiamtes } 
#'   \item{Question}{IPCC question number} 
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/19207697/}
"IPCC"

#' IPCC data-set - Wide format
#'
#' The IPCC-wide data-set  comprises the best estimates 
#' for the phrases "likely" and "unlikely" in six IPCC report sentences. 
#'
#' @format A data frame with 4014 rows and 8 variables:
#' \describe{
#'   \item{Q4}{Each column indicates the estimates for one sentence.}
#'   \item{Q5}{Each column indicates the estimates for one sentence.}
#'   \item{Q6}{Each column indicates the estimates for one sentence.}
#'   \item{Q8}{Each column indicates the estimates for one sentence.}
#'   \item{Q9}{Each column indicates the estimates for one sentence.}
#'   \item{Q10}{Each column indicates the estimates for one sentence.}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/19207697/}
"IPCC_Wide"

#' IPCC data-set - Australian data
#'
#' The IPCC-AUS data-set  comprises the best estimates for the phrases 
#' in IPCC report sentences. 
#'
#' @format A data frame with 4014 rows and 8 variables:
#' \describe{
#'   \item{ID}{Subject ID}
#'   \item{gender}{Gender of subjects, `0`is male, `1`is female}
#'   \item{age}{age of subjects}
#'   \item{cfprob}{personal probability.}
#'   \item{bestprob}{nominated probability.}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/19207697/}
"IPCCAUS"


#' Extinction Study data-set
#'
#' Probability of Human Extinction Study
#'
#' @format A data frame with 1170 rows and 11 variables:
#' \describe{
#'   \item{ID}{Subject ID}
#'   \item{gend}{Gender of subjects, `0`is male, `1`is female}
#'   \item{nation}{The nation of the participants come from}
#'   \item{UK}{effect coding for nation}   
#'   \item{IND}{effect coding for nation}  
#'   \item{political}{political orientation of subjects}
#'   \item{format}{The format of probability elicitation}
#'   \item{order}{the order of probability judgement task.}
#'   \item{SECS_6}{Social conservativsm question on attitude toward gun ownership.}   
#'   \item{EQ1_P}{Probability estimates for general threats.}
#'   \item{EQ3_P}{Probability estimates for the greatest threat.}
#' }
#' @source \url{https://www.michaelsmithson.online/}
"ExtEvent"

#' Patient Time Data
#'
#' Data from Modeling Proportion of Patient Time in Emergency Ward Stages
#'
#' @format A data frame with 1170 rows and 11 variables:
#' \describe{
#'   \item{id}{case identification}
#'   \item{Day}{day of the week ( 0 = Sunday)}
#'   \item{Ambulance}{0 = walk-in; 1 = ambulance-arrival}
#'   \item{Triage}{triage level}   
#'   \item{Triage1}{1 = triage level 1}  
#'   \item{Triage2}{1 = triage level 2}  
#'   \item{Triage3}{1 = triage level 3}  
#'   \item{Triage4}{1 = triage level 4}  
#'   \item{Triage5}{1 = triage level 5}  
#'   \item{Lab}{1 = laboratory test(s) conducted}  
#'   \item{Xray}{1 = x-ray conducted}
#'   \item{Other}{1 = other intervention}
#'   \item{LOS}{length of stay in minutes}
#'   \item{LOSh}{length of stay in hours}
#'   \item{preg}{proportion of time in registration stage}
#'   \item{ptriage}{proportion of time in triage stage}
#'   \item{pnurse}{proportion of time in nursing care stage}
#'   \item{pphysician}{proportion of time in consultation with physician(s)}
#'   \item{pdecis}{proportion of time in decisional stage}
#'   \item{pregptriage}{preg + ptriage}
#'   \item{pphysdecis}{pphysician + pdecis}
#'   \item{prnurse}{pnurse/(pnurse + pregptriage)}
#'   \item{prphysdec}{pphysdecis /(pphysdecis + pregptriage)}
#' }
#' @source \doi{10.1017/S1481803500006539}
"yoon"

