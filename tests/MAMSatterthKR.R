require(lmerTest)

lm.pred <- lm(as.formula(paste("Lightlevel", "~", 
                               paste(c("TVset","Picture"), collapse="*"), sep="")),
              data=TVbo)
TVbo$x <- scale(predict(lm.pred), scale=FALSE)
lmerTVpic <- lmer(Lightlevel ~ TVset*Picture +   Assessor:x  + (1|Assessor) +
                    (1|TVset:Assessor) + (1|Picture:Assessor) + 
                    (1|TVset:Picture:Assessor), data=TVbo)

## TODO: check with SAS dfs for Satterthwaite and KR to agree
tools::assertWarning(anova(lmerTVpic, type=1)) ## warning: ddf=0 for TVset
anova(lmerTVpic, type=1, ddf="Kenward-Roger")
