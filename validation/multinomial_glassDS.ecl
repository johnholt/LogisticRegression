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
glass_ds := DATASET('~thor::testdata::ml::glass_csv', Glass_Rec, CSV(HEADING(1)));
EXPORT multinomial_glassDS := glass_ds;
