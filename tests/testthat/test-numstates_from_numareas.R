context("Unit-testing numstates_from_numareas in BioGeoBEARS with the testthat package")

#
# Following example testthat usages at e.g.
# http://kbroman.org/pkg_primer/pages/tests.html
#


test_that(desc="Check that cladoRcpp version number is >= 0.15", code={

version_number = packageVersion("cladoRcpp")
TF = version_number >= 0.15

if (TF == FALSE)
	{
	txt = 'STOP ERROR inside test_that(desc="Check that cladoRcpp version number is >= 0.15"): the BioGeoBEARS "testthat" tests, located in BioGeoBEARS/tests, require that cladoRcpp have version 0.15 or higher to work. To get the new version, try "devtools::install_github(repo="nmatzke/cladoRcpp", quick=TRUE, dependencies=FALSE, build_vignettes=FALSE, keep_source=TRUE, local=FALSE, force=TRUE)".'
	
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	}

expect_equal(object=TF, expected=TRUE)

}) # END test_that




test_that(desc="Check that phytools is installed", code={

TF = is.element("phytools", installed.packages()[,1])

if (TF == FALSE)
	{
	txt = 'STOP ERROR inside test_that(desc="Check that phytools is installed"): the BioGeoBEARS "testthat" tests, located in BioGeoBEARS/tests, require that phytools be installed. To get it, try \n\ninstall.packages("devtools")\n.\n'
	
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	}

expect_equal(object=TF, expected=TRUE)

}) # END test_that





test_that(desc="check numstates_from_numareas", code={

numstates = numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
expect_equal(object=numstates, expected=16)

numstates = numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
expect_equal(object=numstates, expected=15)

numstates = numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
expect_equal(object=numstates, expected=15)

numstates = numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)
expect_equal(object=numstates, expected=11)

numstates = numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)
expect_equal(object=numstates, expected=1024)

numstates = numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
expect_equal(object=numstates, expected=56)

}) # END test_that









test_that(desc="check numstates_from_numareas for 16 states", code={

numstates = numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
expect_equal(object=numstates, expected=16)

}) # END test_that

test_that(desc="check numstates_from_numareas for 1024 states", code={

numstates = numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)
expect_equal(object=numstates, expected=1024)

}) # END test_that

test_that(desc="check numstates_from_numareas for 1831 states", code={

numstates = numstates_from_numareas(numareas=20, maxareas=2, include_null_range=TRUE)
expect_equal(object=numstates, expected=1831)

}) # END test_that

