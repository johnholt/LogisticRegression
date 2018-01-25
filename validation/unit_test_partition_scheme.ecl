IMPORT $.^ AS LR;
Types := LR.Types;
Constants := LR.Constants;
Types.Field_Desc gen_fd(UNSIGNED c, UNSIGNED4 card) := TRANSFORM
  SELF.number := (c+1) DIV 2;
  SELF.cardinality := IF(c=1, card, 10);
  SELF := [];
END;
Types.Data_Info gen_di(UNSIGNED c) := TRANSFORM
  card_sw := ((c-1) DIV 6) + 1;
  dim_sw := ((c-1) % 6) + 1;
  card := CHOOSE(card_sw, 10, 20, 40, 100, 500, 1000, 10000);
  SELF.wi := c;
  SELF.independent_records := (((c-1) DIV 6)+1)*1000000;
  SELF.dependent_records := (((c-1) DIV 6)+1)*1000000;
  SELF.independent_fields := CHOOSE(dim_sw, 9, 39, 99, 399, 799, 9999);
  SELF.dependent_stats := DATASET(3, gen_fd(COUNTER, card));
  SELF := [];
END;
ds := DATASET(7*6, gen_di(COUNTER));

EXPORT unit_test_partition_scheme := LR.Softmax.Partition_Scheme(ds);
