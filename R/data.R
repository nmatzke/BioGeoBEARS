#######################################################
# Data objects that can be called with the data() command
# Google: how to set up data() in R package
#######################################################


setup_of_a_data_object='
# 1. Create or load your data frame / object in the R console
Psychotria_ML_DEC = resDEC

# 2. Automatically save it as an .rda file in the data/ folder
wd = "~/GitHub/BioGeoBEARS/data/"
setwd(wd)
usethis::use_data(Psychotria_ML_DEC, overwrite=TRUE)
list.files()

# (documentation, add apostrophe:
# Title of Your Dataset
#
# A brief description of what the dataset contains and its purpose.
#
# @format A data frame with 5 rows and 2 variables:
# \\describe{
#   \\item{x}{Integer values representing a sequence.}
#   \\item{y}{Character values representing letters.}
# }
# @source \\url{https://yourdatasource.com}
# @examples
# data(my_dataset)
"my_dataset"

# )
' # END setup_of_a_data_object='


#' Psychotria ML DEC results
#'
#' Results of Maximum Likelihood inference under the 
#' BioGeoBEARS DEC inference
#'
#' @format BioGeoBEARS_results_object
#' \describe{
#'   \item{computed_likelihoods_at_each_node}{Summed likelihoods at each node. Calculation: rowSums(Psychotria_ML_DEC$condlikes_of_each_state)}
#'   \item{relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}{Downpass state probabilities normalized to probabilities (branch tops)}
#'   \item{condlikes_of_each_state}{Downpass likelihoods of each state at each node (branch tops)}
#'   \item{relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS}{Downpass state likelihoods normalized to probabilities (branch bottoms)}
#'   \item{relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS}{Uppass state likelihoods normalized to probabilities (branch bottoms)}
#'   \item{relative_probs_of_each_state_at_branch_top_AT_node_UPPASS}{Uppass state probabilities (branch tops)}
#'   \item{ML_marginal_prob_each_state_at_branch_bottom_below_node}{Marginal probabilities of each state (combining downpass and uppass), at branch bottoms}
#'   \item{ML_marginal_prob_each_state_at_branch_top_AT_node}{Marginal probabilities of each state (combining downpass and uppass), at nodes/branch tops}
#'   \item{relative_probs_of_each_state_at_bottom_of_root_branch}{Downpass state likelihoods normalized to probabilities (bottom of root branch; very rarely used)}
#'   \item{total_loglikelihood}{Total log likelihood. Calculation: sum(log(Psychotria_ML_DEC$computed_likelihoods_at_each_node))}
#'   \item{inputs}{The BioGeoBEARS_run_object input into bears_optim_run() for the ML search}
#'   \item{outputs}{BioGeoBEARS_model_object with the params_table with the updated "est" column for the estimated parameters. Access with Psychotria_ML_DEC$outputs@params_table}
#'   \item{optim_result}{The result from optim, optimx, or GenSA search}
#' }
#' @source \url{https://phylo.wikidot.com/biogeobears#script}
#' @examples
#' data(Psychotria_ML_DEC)
"Psychotria_ML_DEC"
