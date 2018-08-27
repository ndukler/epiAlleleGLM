library(epiAllele)
## Setup dataset for testing
## Create tree
tree=ape::read.tree(text = "((A,B),C);")
tree=ape::unroot.phylo(tree)
tree=ape::reorder.phylo(tree,"postorder")
tree$edge.length=c(0.25,0.5,2)

## Settings
nAlleles=2
species=c("A","B","C")
rate=c(0.25)
pi=c(0.25,0.75)
states=matrix(c(0,1,1),nrow=1)
colnames(states)=tree$tip.label

## Compute site probabilities using epiAllele and construct alleleData object
aData=disCharToProb(states,c(0,1))
ad=alleleData(data=aData,tree=tree)

## Construct two edge group edgeTable
et=getEdgeTable(ad)
et[,edgeGroup:=c("e1","e2","e3")]

## construct rateModel 
rateMod=rateModel(ad,rate = rate,pi = pi,lineageTable = et)

## Set rates on each lineage
epiAllele:::setParamValue(rateMod,getRateIndex(rateMod,edgeGroup = "e1",siteLabel = "All"),value = 1)
epiAllele:::setParamValue(rateMod,getRateIndex(rateMod,edgeGroup = "e2",siteLabel = "All"),value = 0.8)
epiAllele:::setParamValue(rateMod,getRateIndex(rateMod,edgeGroup = "e3",siteLabel = "All"),value = 0.5)

##### Begin Tests ######
## -------------------------------------------------------------------------- ##
testthat::context("rateModel object getter/setter functions")
## Check that values get set correctly
testthat::test_that("Correct rates are retrieved using edgeGroup",
                    testthat::expect_equal(getParamValue(obj = rateMod,getRateIndex(rateMod,edgeGroup = "e1",siteLabel = "All")),1))
testthat::test_that("Correct rates are retrieved using edge table",
                    testthat::expect_equal(getParamValue(obj = rateMod,getRateIndex(rateMod,edges = data.frame(parent=4,child=3),
                                                                                    siteLabel = "All")),0.5))
## Check that the stationary distribution of nucleotide frequencies are retrieved correctly
testthat::test_that("Correct stationary distibution is retrieved",
                    testthat::expect_equal(getParamValue(obj = rateMod,getPiIndex(rateMod,siteLabel = "All")),c(0.25,0.75)))

## -------------------------------------------------------------------------- ##
testthat::context("rateModel object log-likelihood calculations")

## Compute transition tables by hand and compare to function computed ones
qBase=matrix(c(-0.75,0.75,0.25,-0.25),ncol=2,nrow = 2,byrow = TRUE)
qBaseNorm=qBase/(2*0.75*0.25)
qBaseNormE1=log(as.matrix(Matrix::expm(qBaseNorm*0.25*1)))
qBaseNormE2=log(as.matrix(Matrix::expm(qBaseNorm*0.5*0.8)))
qBaseNormE3=log(as.matrix(Matrix::expm(qBaseNorm*2*0.5)))

## Compute tables using functions
tt=data.table::data.table(getTree(rateMod)$edge)
data.table::setnames(x = tt,old=colnames(tt),c("parent","child"))
rates=getParams(rateMod)[getRateIndex(rateMod,edges = tt,siteLabel = "All")]
pi=getParams(rateMod)[getPiIndex(rateMod,siteLabel = "All")]
## Compute transition matricies
logTransMat=epiAllele:::branchRateMatrix(rate = rates,branch.length =  getTree(rateMod)$edge.length,pi = pi)
## Compute log-likelihood
pl1=exp(qBaseNormE1) %*% matrix(c(1,0),ncol = 1)
pl2=exp(qBaseNormE2) %*% matrix(c(0,1),ncol=1)  
pl3=exp(qBaseNormE3) %*% matrix(c(0,1),ncol=1)  

## Check that the stationary distribution of nucleotide frequencies are retrieved correctly
testthat::test_that("The correct log transition matrix is computed - E1",
                    testthat::expect_equal(logTransMat[[1]],qBaseNormE1))
testthat::test_that("The correct log transition matrix is computed - E2",
                    testthat::expect_equal(logTransMat[[2]],qBaseNormE2))
testthat::test_that("The correct log transition matrix is computed - E3",
                    testthat::expect_equal(logTransMat[[3]],qBaseNormE3))

## Check that the manually computed log-likelihood is equal to the output of the function
testthat::test_that("The correct log-likelihood is computed - E3",
                    testthat::expect_equal(logSumExp(log(pl1)+log(pl2)+log(pl3)+log(pi)),logLikelihood(obj=rateMod)))

