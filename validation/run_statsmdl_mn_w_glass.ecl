// Run python Stats Model Multinomial Logistic Regression
IMPORT Python;
Glass_Rec := RECORD
  REAL8 RI;
  REAL8 Na;
  REAL8 Mg;
  REAL8 Al;
  REAL8 Si;
  REAL8 K;
  REAL8 Ca;
  REAL8 Ba;
  REAL8 Fe;
  INTEGER Type_Glass;
END;
glass_ds := DATASET('~thor::jdh::test::glass.csv', Glass_Rec, CSV(HEADING(1)));
//t0 := TABLE(glass_ds, {Type_Glass, c:=COUNT(GROUP)}, Type_Glass, FEW, UNSORTED);
//UNSIGNED statsmdl_mn(DATASET(Glass_Rec) g) := EMBED(Python)
UNSIGNED statsmdl_mn(UNSIGNED g) := EMBED(Python)
count(g)
ENDEMBED;
//statsmdl_mn(glass_ds);
statsmdl_mn(7);
//EXPORT run_statsmdl_mn_w_glass := OUTPUT(t0, NAMED('Python_rslt'));
