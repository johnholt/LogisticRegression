IMPORT $.^ AS LR;
IMPORT LR.Constants;
IMPORT LR.Types;
IMPORT $ AS Smax;
IMPORT Std.System.ThorLib;
// aliases for convenience
Info := Types.Data_Info;
Scheme := Types.Part_Scheme;
t_work_item := Types.t_work_item;
NumericField := Types.NumericField;
/**
 * Determine the set of partitioning schemes for the collection of work
 * items.  The Constants.local_cap value is used as a ceiling.
 * @param ds The data information dataset.
 * @returns partitioning scheme for each work item
 */
EXPORT DATASET(Scheme) Partition_Scheme(DATASET(Info) ps) := FUNCTION
  ps2 := ASSERT(ps, EXISTS(dependent_stats(number=1 AND cardinality>0)),
                'Dependents not col 1 or no card, run Data_Check', FAIL);
  Scheme ext1(Info i, Types.Field_Desc f) := TRANSFORM
    obs := MAX(i.dependent_records, i.independent_records);
    dim := i.independent_fields + 1;  // plus 1 for intercept
    cls := f.cardinality;
    cols:= MAX(dim, cls);
    // Figure class-dims matrix first
    coef_cells := dim * cls;
    coef_parts := ((coef_cells-1) DIV Constants.local_cap) + 1;
    dim_part := ((coef_parts*dim) DIV (dim+cls)) + 1;
    cls_part := MAX(1, coef_parts-dim_part);
    //  Now figure data partitioning
    col_part := MAX(cls_part, dim_part);
    data_cells := obs*cols;
    data_parts := ((data_cells-1) DIV (Constants.local_cap*col_part)) + 1;
    //
    SELF.wi := i.wi;
    SELF.obs := obs;
    SELF.cls := cls;
    SELF.dim := dim;
    SELF.obs_per_part := ((obs-1) DIV data_parts) + 1;
    SELF.dim_per_part := ((dim-1) DIV dim_part) + 1;
    SELF.cls_per_part := ((cls-1) DIV cls_part) + 1;
    SELF.data_nodes := MIN(((data_parts-1) DIV col_part)+1, ThorLib.nodes());
    SELF.data_parts := data_parts;
    SELF.coef_parts := coef_parts;
    SELF.coef_nodes := MIN(coef_parts, ThorLib.nodes());;
    SELF.frst_node := 0;
  END;
  pick1(DATASET(Types.Field_Desc) ds) := CHOOSEN(ds(number=1), 1);
  t0 := NORMALIZE(ps2, pick1(LEFT.dependent_stats), ext1(LEFT, RIGHT));
  Scheme enum_nodes(Scheme prev, Scheme curr) := TRANSFORM
    SELF.frst_node := (prev.frst_node + prev.data_nodes) % ThorLib.nodes();
    SELF := curr;
  END;
  t1 := ITERATE(t0, enum_nodes(LEFT, RIGHT));
  RETURN t1;
END;
