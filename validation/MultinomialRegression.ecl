// Multinomial Regression validation
IMPORT $ AS V;
IMPORT $.^ AS LR;
IMPORT ML_Core;
IMPORT ML_Core.Types AS Core_Types;
IMPORT LR.Types AS Types;

TestData := RECORD
  Types.t_recordID id;
  RECORDOF(V.multinomial_glassDS);
END;
testData fix_recs(V.multinomial_glassDS rec, UNSIGNED c) := TRANSFORM
  SELF.id := c;
  SELF.type_glass := IF(rec.type_glass > 3, -1, 0) + rec.type_glass;
  SELF := rec;
END;
glass := PROJECT(V.multinomial_glassDS, fix_recs(LEFT,COUNTER));
glass_even := PROJECT(glass(id%2=0), TRANSFORM(TestData, SELF.id:=COUNTER,SELF:=LEFT));
glass_odd := PROJECT(glass(id%2=1), TRANSFORM(TestData, SELF.id:=COUNTER,SELF:=LEFT));

ML_Core.ToField(glass, glass_indep_all, id, , 1,
             'RI,Na,Mg,Al,Si,K,Ca,Ba,Fe');
ML_Core.ToField(glass_even, glass_indep_even, id, , 2,
             'RI,Na,Mg,Al,Si,K,Ca,Ba,Fe');
ML_Core.ToField(glass_odd, glass_indep_odd, id, , 3,
             'RI,Na,Mg,Al,Si,K,Ca,Ba,Fe');
glass_indep := glass_indep_all + glass_indep_even + glass_indep_odd;
ML_Core.ToField(glass, glass_dep_all, id, , 1, 'Type_Glass');
ML_Core.ToField(glass_even, glass_dep_even, id, , 2, 'Type_Glass');
ML_Core.ToField(glass_odd, glass_dep_odd, id, , 3, 'Type_Glass');
glass_classes := PROJECT(glass_dep_all+glass_dep_even+glass_dep_odd, Types.DiscreteField);

mod := LR.Softmax.GetModel(glass_indep, glass_classes, 2000);

EXPORT MultinomialRegression := OUTPUT(mod(wi=1), ALL);