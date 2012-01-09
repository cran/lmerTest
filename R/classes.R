merLmerTest<-setClass("merLmerTest", representation(t.pval="numeric") ,contains = "mer")
summary.merLmerTest<-setClass("summary.merLmerTest", contains = "summary.mer")