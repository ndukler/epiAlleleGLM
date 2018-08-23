context("treeLL")

## Create tree
tree=ape::read.tree(text = "((A,B),C);")
tree=ape::unroot.phylo(tree)
tree=ape::reorder.phylo(tree,"postorder")
tree$edge.length=c(0.25,0.5,1)

## Settings
nAlleles=2
species=c("A","B","C")
rate=c(0.25)
subsetStates=1

## Compute stationary dist probabilities
pi=c(0.25,0.75)

## Enumerate all states and get relevant subset
states=expand.grid(lapply(as.list(1:length(tree$tip.label)),function(x) c(1:nAlleles)))[1,]
colnames(states)=tree$tip.label

## Compute site probabilities using epiAllele and construct alleleData object
aData=lapply(split(t(states), f =colnames(states)),function(x){
  z=matrix(0,nrow = length(x),ncol=nAlleles)
  for(k in 1:nrow(z)){z[k,x[k]]=1}
  return(z)
})
ad=alleleData(data=aData,tree=tree)

## Construct two edge group edgeTable
et=getEdgeTable(ad)
et[,edgeGroup:=c("e1","e1","e2")]

## construct rateModel 
rateMod=rateModel(ad,rate = rate,pi = pi,lineageTable = et)

epiAllele:::setParamValue(rateMod,getRateIndex(rateMod,edgeGroup = "e2",siteLabel = "All"),value = 2)
plotTree(rateMod,"value")

eLL=logLikelihood(obj=rateMod)







## Compute site probabilities using phangorn
statesP=phangorn::phyDat(states,type="USER",levels=c(1:nAlleles))
pLL=phangorn::pml(tree,data=statesP,bf=pi,rate=rate)

## Compute site probabilities using epiAllele
aData=lapply(split(t(states), f =colnames(states)),function(x){
  z=matrix(0,nrow = length(x),ncol=nAlleles)
  for(k in 1:nrow(z)){z[k,x[k]]=1}
  return(z)
})

ad=alleleData(data=aData,tree=tree)
rateMod=rateModel(ad,rate = rate,pi = pi)
eLL=logLikelihood(obj=rateMod)


pLL$logLik
eLL
