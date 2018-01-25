IMPORT ML_Core;
IMPORT ML_Core.Types AS Core_Types;
IMPORT $.^ AS LR;
IMPORT LR.Constants;
IMPORT LR.Types;
IMPORT $ AS Smax;
IMPORT Std.System.ThorLib;
IMPORT Std.BLAS;
//Aliases for convenience
AnyField := Types.AnyField;
NumericField := Types.NumericField;
DiscreteField := Types.DiscreteField;
Layout_Model := Types.Layout_Model;
t_work_item  := Types.t_work_item;
t_RecordID  := Types.t_RecordID;
t_FieldNumber := Types.t_FieldNumber;
t_FieldReal := Types.t_FieldReal;
t_Discrete := Types.t_Discrete;
t_dimension := Types.t_dimension;
t_matrix := Types.t_matrix;
Data_Info := Types.Data_Info;
Part_Scheme := Types.Part_Scheme;
Part_Rec := Types.Part_Rec;
t_value      := BLAS.Types.value_t;
triangle     := BLAS.Types.triangle;
diagonal     := BLAS.Types.diagonal;
side         := BLAS.Types.side;
Apply2Cells  := BLAS.Apply2Cells;
gemm         := BLAS.dgemm;
trsm         := BLAS.dtrsm;
getf2        := BLAS.dgetf2;
potf2        := BLAS.dpotf2;
axpy         := BLAS.daxpy;
asum         := BLAS.dasum;
scal         := BLAS.dscal;
make_diag    := BLAS.make_diag;
make_vector  := BLAS.make_vector;
extract_diag := BLAS.extract_diag;
dimm         := LR.dimm;
gemagg       := LR.gemagg;

/**
 * Generate logistic regression model from training data.  The size
 * of the inputs is used to determine which work items are processed
 * with purely local operations (the data is moved once as necessary)
 * or with global operations supporting a work item to use multiple
 * nodes.
 * The classes by dimensions matrix of co-efficients must fit in a single
 * partition.
 * The training data classe names must be dense and start at 1.
 * @param independents the independent values
 * @param dependents the dependent values.  There should be only 1 dependent.
 * @param max_iter maximum number of iterations to try
 * @param epsilon the minimum change in the Theta value estimate to continue
 * @param alpha the learning rate, default is 0.1
 * @param lambda the decay rate, default is 0.0001
 * @return coefficient matrix plus model building stats
 */
EXPORT //DATASET(Layout_Model)
      GetModel(DATASET(NumericField) independents,
               DATASET(DiscreteField) dependents,
               UNSIGNED max_iter=200,
               REAL8 epsilon=Constants.default_epsilon,
               REAL8 alpha=0.1, REAL8 lambda=0.0001) := FUNCTION
  //**************************************************
  //**************************************************
  // Pick up data statistics
  //**************************************************
  stats := LR.DataStats(independents, dependents, FALSE, TRUE);
  s0 := Smax.Partition_Scheme(stats);
  scheme := s0(dim=dim_per_part AND cls=cls_per_part);
  initial_theta := Smax.InitialTheta(scheme);
  //**************************************************
  //**************************************************
  // Conversion to dense partitions from sparse cells
  //**************************************************
  NF_Work := RECORD(NumericField)
    t_dimension matrix_rows;
    t_dimension matrix_cols;
    t_dimension part_rows;
    t_dimension part_cols;
    t_dimension node_id;
    t_dimension block_row;
    t_dimension block_col;
    t_dimension first_node_id;
    t_dimension part; // zero based
    t_dimension parts;
    UNSIGNED4 inserts;
    BOOLEAN dropMe;
  END;
  NF_Work gen_dummy(Part_Scheme s, UNSIGNED c, t_dimension cols) := TRANSFORM
    part := (c-1) % s.data_parts;
    SELF.id := MIN((part+1)*s.obs_per_part, s.obs);
    SELF.number := ((c-1) DIV s.data_parts) + 1;
    SELF.value := 0;
    SELF.wi := s.wi;
    SELF.matrix_rows := s.obs;
    SELF.matrix_cols := cols;
    SELF.part_rows := MIN(s.obs_per_part, s.obs-(part*s.obs_per_part));
    SELF.part_cols := cols;
    SELF.node_id := (s.frst_node + part) % ThorLib.nodes();
    SELF.block_row := part + 1;
    SELF.block_col := 1;
    SELF.first_node_id := s.frst_node;
    SELF.part := part;
    SELF.parts := s.data_parts;
    SELF.inserts := 0;
    SELF.dropMe := TRUE;
  END;
  NF_Work set_inserts(NF_Work prev, NF_Work curr) := TRANSFORM
    SELF.inserts := curr.id - prev.id -1;
    SELF.dropMe := curr.dropMe AND curr.id=prev.id;
    SELF := curr;
  END;
  NF_Work make_inserts(NF_Work nf, UNSIGNED c) := TRANSFORM
    SELF.value := IF(c>nf.inserts, nf.value, 0.0);
    SELF.id := nf.id - nf.inserts + c - 1;
    SELF := nf;
  END;
  Value_Rec := {t_FieldReal value};
  Part_Rec roll2part(NF_Work nf, DATASET(NF_Work) nfs, BOOLEAN insert) := TRANSFORM
    ones := DATASET(nf.part_rows, TRANSFORM(Value_Rec, SELF.value:=1.0));
    mat_ds := IF(insert, ones)  & PROJECT(nfs, Value_Rec);
    SELF.mat := SET(mat_ds, value);
    SELF := nf;
  END;
  //**************************************************
  //**************************************************
  // Convert indep NumericFormat to partition records, adding 1 column
  //**************************************************
  NF_Work cvt_indep(NumericField nf, Part_Scheme s) := TRANSFORM
    part := (nf.id-1) DIV s.obs_per_part;
    SELF.matrix_rows := s.obs;
    SELF.matrix_cols := s.dim;
    SELF.part_rows := MIN(s.obs_per_part, s.obs-(part*s.obs_per_part));
    SELF.part_cols := s.dim;
    SELF.node_id := (s.frst_node + part) % ThorLib.nodes();
    SELF.block_row := part + 1;
    SELF.block_col := 1;
    SELF.first_node_id := s.frst_node;
    SELF.part := part;
    SELF.parts := s.data_parts;
    SELF.inserts := 0;
    SELF.dropMe := FALSE;
    SELF := nf;
  END;
  scheme_ind := JOIN(independents, scheme, LEFT.wi=RIGHT.wi,
                     cvt_indep(LEFT, RIGHT), LOOKUP);
  dummy_ind := NORMALIZE(scheme, LEFT.data_parts*LEFT.dim,
                         gen_dummy(LEFT, COUNTER, LEFT.dim));
  dist_ind := DISTRIBUTE(scheme_ind+dummy_ind, node_id);
  srtd_ind := SORT(dist_ind, wi, part, number, id, dropMe, LOCAL);
  grpd_ind := GROUP(srtd_ind, wi, part, number, LOCAL);
  marked_ind := ITERATE(grpd_ind, set_inserts(LEFT, RIGHT));
  used_ind := marked_ind(NOT dropMe);
  dense_ind := NORMALIZE(used_ind, LEFT.inserts+1, make_inserts(LEFT, COUNTER));
  ready_ind := GROUP(dense_ind, wi, part, LOCAL);
  mat_ind := ROLLUP(ready_ind, GROUP, roll2part(LEFT, ROWS(LEFT), TRUE));
  //**************************************************
  //**************************************************
  // Convert dep Discrete Format to partition records as ground truth
  //**************************************************
  NF_Work cvt_dep(DiscreteField df, Part_Scheme s) := TRANSFORM
    part := (df.id-1) DIV s.obs_per_part;
    SELF.matrix_rows := s.obs;
    SELF.matrix_cols := s.cls;
    SELF.part_rows := MIN(s.obs_per_part, s.obs-(part*s.obs_per_part));
    SELF.part_cols := s.cls;
    SELF.node_id := (s.frst_node + part) % ThorLib.nodes();
    SELF.block_row := part + 1;
    SELF.block_col := 1;
    SELF.first_node_id := s.frst_node;
    SELF.part := part;
    SELF.parts := s.data_parts;
    SELF.inserts := 0;
    SELF.dropMe := FALSE;
    SELF.number := df.value;
    SELF.value := 1.0;
    SELF := df;
  END;
  scheme_dep := JOIN(dependents(value>0), scheme, LEFT.wi=RIGHT.wi,
                     cvt_dep(LEFT, RIGHT), LOOKUP);
  dummy_dep := NORMALIZE(scheme, LEFT.data_parts*LEFT.cls,
                         gen_dummy(LEFT, COUNTER, LEFT.cls));
  dist_dep := DISTRIBUTE(scheme_dep+dummy_dep, node_id);
  srtd_dep := SORT(dist_dep, wi, part, number, id, dropMe, LOCAL);
  grpd_dep := GROUP(srtd_dep, wi, part, number, LOCAL);
  marked_dep := ITERATE(grpd_dep, set_inserts(LEFT, RIGHT));
  used_dep := marked_dep(NOT dropMe);
  dense_dep := NORMALIZE(used_dep, LEFT.inserts+1, make_inserts(LEFT, COUNTER));
  ready_dep := GROUP(dense_dep, wi, part, LOCAL);
  mat_dep := ROLLUP(ready_dep, GROUP, roll2part(LEFT, ROWS(LEFT), FALSE));
  //**************************************************
  //**************************************************
  // Iteration function definitions
  //**************************************************
  t_value e(t_value v, t_dimension r, t_dimension c) := EXP(v);
  t_value reciprocal(t_value v, t_dimension r, t_dimension c) := 1.0/v;
  t_value absval(t_value v, t_dimension r, t_dimension c) := ABS(v);
  Part_Rec class_probs(Part_Rec x_part, Part_Rec t_part) := TRANSFORM
    cls := t_part.matrix_cols;
    dim := t_part.matrix_rows;
    part_obs := x_part.part_rows;
    ones := make_vector(cls, 1.0);
    obs_score := gemm(FALSE, FALSE, part_obs, cls, dim,
                      1.0, x_part.mat, t_part.mat);
    max_score := gemagg(part_obs, cls, obs_score,
                        FALSE, Types.Matrix_Op.Max_Val);
    adj_score := gemm(FALSE, FALSE, part_obs, cls, 1,
                      -1.0, max_score, ones, 1.0, obs_score);
    exp_adj := Apply2Cells(part_obs, cls, adj_score, e);
    sum_exp_adj := gemagg(part_obs, cls, exp_adj,
                          FALSE, Types.Matrix_Op.Sum_Val);
    recip_sum := Apply2Cells(part_obs, 1, sum_exp_adj, reciprocal);
    SELF.matrix_cols := cls;
    SELF.part_cols := cls;
    SELF.mat := dimm(FALSE, FALSE, TRUE, FALSE, part_obs, cls, part_obs,
                     1.0, recip_sum, exp_adj);
    SELF := x_part;
  END;
  Part_Rec calc_term(DATASET(Part_Rec) parts) := TRANSFORM
    y_part := parts[1];
    p_part := parts[2];
    x_part := parts[3];
    part_obs := y_part.part_rows;
    cls := y_part.part_cols;
    dim := x_part.part_cols;
    G_less_P := axpy(part_obs*cls, -1.0, p_part.mat, 1, y_part.mat, 1);
    SELF.mat := gemm(TRUE, FALSE, dim, cls, part_obs,
                     1.0, x_part.mat, G_less_P);
    SELF.matrix_rows := dim; // Dim x Cls
    SELF.matrix_cols := cls;
    SELF.part_rows := dim;
    SELF.part_cols := cls;
    SELF.block_row := 1;
    SELF.block_col := 1;
    SELF.parts := 1;
    SELF.part := 0;
    SELF.node_id := y_part.first_node_id;
    SELF := y_part;
  END;
  Part_Rec sum_parts(Part_Rec base, Part_Rec incr) := TRANSFORM
    dim := base.matrix_rows;
    cls := base.matrix_cols;
    SELF.mat := axpy(dim*cls, 1.0, base.mat, 1, incr.mat, 1);
    SELF := base;
  END;
  Iter_Part := RECORD(Part_Rec)
    UNSIGNED4 iter;
    REAL8 delta;
  END;
  Iter_Part add_fields(Part_Rec p) := TRANSFORM
    SELF.iter := 0;
    SELF.delta := 2*epsilon;
    SELF := p;
  END;
  Iter_Part update_theta(DATASET(Iter_Part) rws) := TRANSFORM
    upd := rws[1];
    old := rws[2];
    dat := rws[3];  // needed to get number of observations
    obs := dat.matrix_rows;
    cls := upd.matrix_cols;
    dim := upd.matrix_rows;
    decayed := scal(dim*cls, (1.0-alpha*lambda), old.mat, 1);
    new_mat := axpy(dim*cls, (alpha/obs), upd.mat, 1, decayed, 1);
    diff := axpy(dim*cls, -1.0, old.mat, 1, new_mat, 1);
    abs_diff := Apply2Cells(dim, cls, diff, absval);
    SELF.iter := old.iter + 1;
    SELF.delta := MAX(abs_diff);
    SELF.mat := new_mat;
    SELF := old;
  END;
  DATASET(Iter_Part) bfgs(DATASET(Iter_Part) T) := FUNCTION
    Cls_Probs_Est := JOIN(mat_ind, T, LEFT.wi=RIGHT.wi,
                          class_probs(LEFT, RIGHT), LOOKUP);
    data_parts := [mat_dep, Cls_Probs_Est, mat_ind];
    partial_terms := JOIN(data_parts,
                          LEFT.wi=RIGHT.wi AND LEFT.part=RIGHT.part,
                          calc_term(ROWS(LEFT)),
                          SORTED(wi, part), LOCAL);
    dist_parts := DISTRIBUTE(partial_terms, first_node_id);
    ready_parts := SORT(dist_parts, wi, part, LOCAL);
    iter_adj := ROLLUP(ready_parts, sum_parts(LEFT, RIGHT), wi, LOCAL);
    theta_parts := [PROJECT(iter_adj, add_fields(LEFT)),
                    T, PROJECT(mat_ind(part=0), add_fields(LEFT))];
    newT := JOIN(theta_parts, LEFT.wi=RIGHT.wi,
                 update_theta(ROWS(LEFT)), SORTED(wi), LOCAL);
    RETURN newT;
  END;
  //**************************************************
  //**************************************************
  // Loop until converged or reached max iterations
  //**************************************************
  theta := LOOP(PROJECT(initial_theta, add_fields(LEFT)), max_iter,
                LEFT.iter<max_iter AND ABS(LEFT.delta) > epsilon,
                bfgs(ROWS(LEFT)), UNORDERED);
  //**************************************************
  //**************************************************
  // Calculate class probability estimates for SE
  //**************************************************
  All_Cls_Probs := JOIN(mat_ind, theta, LEFT.wi=RIGHT.wi,
                        class_probs(LEFT, RIGHT), LOOKUP);
  Cls_Prob_Part := RECORD(Part_Rec)
    t_FieldNumber this_cls;
  END;
  Cls_Prob_Part extract_cls(Part_Rec pr, UNSIGNED col) := TRANSFORM
    first_pos := ((col-1)*pr.part_rows) + 1;
    last_pos := first_pos + pr.part_rows-1;
    SELF.mat := pr.mat[first_pos..last_pos];
    SELF.this_cls := col;
    SELF.part_cols := 1;
    SELF.matrix_cols := 1;
    SELF := pr;
  END;
  Each_Cls_Prob := NORMALIZE(All_Cls_Probs,LEFT.part_cols, extract_cls(LEFT,COUNTER));
  //**************************************************
  //**************************************************
  // Calculate SE of estimate
  //**************************************************
  Cov_Part := RECORD
    t_work_item wi;
    t_dimension dim;
    t_Discrete this_cls;
    t_matrix mat;    // dim x dim
  END;
  t_value p_not_p(t_value v, t_dimension r, t_dimension c) := v*(1-v);
  Cov_Part calc_cov_term(Part_Rec x, Cls_Prob_Part p) := TRANSFORM
    V_diag := Apply2Cells(p.part_rows, 1, p.mat, p_not_p);
    VX := dimm(FALSE, FALSE, TRUE, FALSE, x.part_rows, x.part_cols, x.part_rows,
               1.0, V_diag, x.mat);
    SELF.mat := gemm(TRUE, FALSE, x.part_cols, x.part_cols, x.part_rows,
                     1.0, x.mat, VX);
    SELF.dim := x.part_cols;
    SELF.this_cls := p.this_cls;
    SELF.wi := x.wi;
  END;
  cov_0 := JOIN(mat_ind, Each_Cls_Prob,
                   LEFT.wi=RIGHT.wi AND LEFT.part=RIGHT.part,
                   calc_cov_term(LEFT,RIGHT), LOCAL);
  cov_terms := SORT(cov_0, wi, this_cls);
  Cov_Part sum_cov_trms(Cov_Part base, Cov_Part incr) := TRANSFORM
    entries := base.dim*base.dim;
    SELF.mat := axpy(entries, 1.0, base.mat, 1, incr.mat, 1);
    SELF := base;
  END;
  cov_mat := ROLLUP(cov_terms, sum_cov_trms(LEFT, RIGHT), wi, this_cls, LOCAL);
  Cov_Part calc_se(Cov_Part cov) := TRANSFORM
    // Solve Ax = I, but that might not be stable
    // so A'Ax = A' is stable and A'A is positive definite if A is full rank
    // Don't have a solve of a transform, so x'A'A = A instead
    // as we only need diagonal entries for SE
    AtA := gemm(TRUE, FALSE, cov.dim, cov.dim, cov.dim, 1.0, cov.mat, cov.mat);
//    AtA_fac := potf2(triangle.upper, cov.dim, AtA, TRUE);
//    AtA_fac2:= potf2(triangle.lower, cov.dim, AtA, TRUE);
    AtA_lu := getf2(cov.dim, cov.dim, AtA);
    intermed := trsm(side.xA, triangle.upper, FALSE, diagonal.NotUnitTri,
                     cov.dim, cov.dim, cov.dim, 1.0, AtA_lu, cov.mat);
    inv_covT := trsm(side.xA, triangle.lower, FALSE, diagonal.UnitTri,
                     cov.dim, cov.dim, cov.dim, 1.0, AtA_lu, intermed);
    SELF.mat := extract_diag(cov.dim, cov.dim, inv_covT);
    //SELF.mat := AtA_fac;
    //SELF.mat := gemm(TRUE, FALSE, cov.dim, cov.dim, cov.dim, 1.0, inv_covT, cov.mat);
    SELF := cov;
  END;
  sse_mat := PROJECT(cov_mat, calc_se(LEFT));
  //**************************************************
  //**************************************************
  // Count score
  //**************************************************
  Score_Part := RECORD
    t_work_item wi;
    UNSIGNED4 correct;
    UNSIGNED4 incorrect;
  END;
  Score_Part score_obs(Part_Rec probs, Part_Rec y) := TRANSFORM
    part_obs := probs.part_rows;
    cls := probs.part_cols;
    max_prob_raw_pos := gemagg(part_obs, cls, probs.mat, FALSE, Types.Matrix_Op.Max_Pos);
    act_class_raw_pos := gemagg(part_obs, cls, y.mat, FALSE, Types.Matrix_Op.Max_Pos);
    diff := axpy(part_obs, -1.0, act_class_raw_pos, 1, max_prob_raw_pos, 1);
    correct := COUNT(DATASET(diff, {REAL8 value})(value=0));
    SELF.correct := correct;
    SELF.incorrect := part_obs - correct;
    SELF.wi := probs.wi;
  END;
  s_part := JOIN(All_Cls_Probs, mat_dep,
                 LEFT.wi=RIGHT.wi AND LEFT.part=RIGHT.part,
                 score_obs(LEFT,RIGHT), LOCAL);
  s_sum := ROLLUP(GROUP(s_part, wi, ALL), GROUP,
                  TRANSFORM(Score_Part, SELF.wi:=LEFT.wi,
                            SELF.correct:=SUM(ROWS(LEFT), correct),
                            SELF.incorrect:=SUM(ROWS(LEFT), incorrect)));
  //**************************************************
  //**************************************************
  // Build model data set
  //**************************************************
  Layout_Model extCoef_SE(Iter_Part p, UNSIGNED subscript) := TRANSFORM
    this_class := ((subscript-1) DIV p.matrix_rows) + 1;
    this_dim := ((subscript-1) % p.matrix_rows) + 1;
    block := Constants.id_betas_coef * p.matrix_rows;
    SELF.wi := p.wi;
    SELF.id := Constants.id_betas + block + this_dim - 1;
    SELF.number := this_class;
    SELF.value := p.mat[subscript];
  END;
  coef := NORMALIZE(theta, LEFT.matrix_rows*LEFT.matrix_cols,
                    extCoef_SE(LEFT, COUNTER));
  Layout_Model extSE(Cov_Part d, UNSIGNED subscript) := TRANSFORM
    block := Constants.id_betas_se*d.dim;
    SELF.wi := d.wi;
    SELF.number := d.this_cls;
    SELF.value := SQRT(d.mat[subscript]);
    SELF.id := Constants.id_betas + block + subscript - 1;
  END;
  se := NORMALIZE(sse_mat, LEFT.dim, extSE(LEFT, COUNTER));
  Layout_Model extScore(Score_Part sp, UNSIGNED c) := TRANSFORM
    SELF.wi := sp.wi;
    SELF.id := IF(c=1, Constants.id_correct, Constants.id_incorrect);
    SELF.number := 1;
    SELF.value := IF(c=1, sp.correct, sp.incorrect);
  END;
  scores := NORMALIZE(s_sum, 2, extScore(LEFT, COUNTER));
  Layout_Model extStats(Iter_Part p, UNSIGNED subscript) := TRANSFORM
    SELF.wi := p.wi;
    SELF.id := IF(subscript=2, Constants.id_delta, Constants.id_iters);
    SELF.number := 1;
    SELF.value := IF(subscript=2, p.delta, p.iter);
  END;
  statistics := NORMALIZE(theta, 2, extStats(LEFT, COUNTER));
  Layout_Model extBase(Data_Info di, UNSIGNED c) := TRANSFORM
    SELF.wi := di.wi;
    SELF.id := Constants.id_base;
    SELF.number := CHOOSE(c, Constants.base_builder,
                             Constants.base_max_iter,
                             Constants.base_epsilon,
                             Constants.base_ind_vars,
                             Constants.base_dep_vars,
                             Constants.base_obs,
                             Constants.base_cls);
    SELF.value := CHOOSE(c, Constants.builder_softmax,
                            max_iter,
                            epsilon,
                            di.independent_fields,
                            COUNT(di.dependent_stats),
                            MAX(di.independent_records, di.dependent_records),
                            di.dependent_stats(number=1)[1].cardinality);
  END;
  base := NORMALIZE(stats, 7, extBase(LEFT, COUNTER));
  rslt := SORT(base+coef+se+scores+statistics, wi, id, number);
  RETURN rslt;
END;