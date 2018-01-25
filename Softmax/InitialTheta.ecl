IMPORT ML_Core;
IMPORT ML_Core.Types AS Core_Types;
IMPORT $.^ AS LR;
IMPORT LR.Constants;
IMPORT LR.Types;
IMPORT $ AS Smax;
IMPORT Std.System.ThorLib;
// aliases for convenience
Part := Types.Part_Rec;
Scheme := Types.Part_Scheme;
t_work_item := Types.t_work_item;
t_dimension := Types.t_dimension;
NumericField := Core_Types.NumericField;
/**
 * Creates a data set of initial co-efficient values.  Initial values
 * are in the range of 0 to 0.005 and are random.  The field number
 * for the value to be predicted is 1.  Dependent records for other fields
 * are ignored.
 * @param s The partition scheme for the data.
 * @returns a Part_Rec matrix of the initial theta values for each work
 * item.  Matrix is Dims x Classes, where dims is the number of dimensions
 * of the independent data plus 1 (for the intercept).
 */
EXPORT DATASET(Part) InitialTheta(DATASET(Scheme) s) := FUNCTION
  NF_ext := RECORD(NumericField)
    t_dimension matrix_rows;
    t_dimension matrix_cols;
    t_dimension first_node_id;
  END;
  NF_ext gen1(Scheme schm, UNSIGNED c) := TRANSFORM
    SELF.wi := schm.wi;
    SELF.id := ((c-1) % schm.dim) + 1;
    SELF.number := ((c-1) DIV schm.dim) + 1;
    SELF.value := 0.005 * (RANDOM()/4294967296.0);
    SELF.matrix_rows := schm.dim;
    SELF.matrix_cols := schm.cls;
    SELF.first_node_id := schm.frst_node;
  END;
  s1 := ASSERT(s, dim_per_part=dim AND cls_per_part=cls,
               'Too many dimensions and classes, run Data_Check', FAIL);
  t0 := SORT(DISTRIBUTE(s1, frst_node), wi, LOCAL);
  c0 := NORMALIZE(t0, LEFT.dim*LEFT.cls, gen1(LEFT, COUNTER));
  // roll up into matrix
  c1 := GROUP(c0, wi, LOCAL);
  Part make_part(NF_ext d, DATASET(NF_ext) rws):=TRANSFORM
    SELF.wi := d.wi;
    SELF.matrix_rows := d.matrix_rows;
    SELF.matrix_cols := d.matrix_cols;
    SELF.part_rows := d.matrix_rows;
    SELF.part_cols := d.matrix_cols;
    SELF.node_id := d.first_node_id;
    SELF.block_row := 1;
    SELF.block_col := 1;
    SELF.first_node_id := d.first_node_id;
    SELF.part := 0;
    SELF.parts := 1;
    SELF.mat := SET(rws, value);
  END;
  rslt := ROLLUP(c1, GROUP, make_part(LEFT, ROWS(LEFT)));
  RETURN rslt;
END;
