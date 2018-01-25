IMPORT $ AS Test;
IMPORT $.^ AS LR;
IMPORT LR.Types;
gemagg := LR.gemagg;
Op := Types.Matrix_Op;

Types.t_matrix t1 := [1,2,3,4,5,6,7,8,0,0,0,1];
Types.t_matrix t2 := [-1,2,3,4,5,6,7,8,0,0,0,1];
M := 4;
N := 3;
Test_In := RECORD
  Types.t_matrix std;
  BOOLEAN use_t1;
  BOOLEAN row_vector;
  Op thisOp
END;
Test_Out := RECORD(Test_In)
  BOOLEAN ok;
  Types.t_matrix tst;
END;
Test_Out ex(Test_In t_rec) := TRANSFORM
  SELF.tst := gemagg(M, N, IF(t_rec.use_t1, t1, t2), t_rec.row_vector, t_rec.thisOp);
  SELF.ok := t_rec.std = SELF.tst;
  SELF := t_rec;
END;
test_set := DATASET([{[4,8,1], TRUE, TRUE, Op.Max_Val}
 ,{[5,6,7,8], TRUE, FALSE, Op.Max_Val}
 ,{[1,5,0], TRUE, TRUE, Op.Min_Val}
 ,{[0,0,0,1], TRUE, FALSE, Op.Min_Val}
 ,{[3,7,11], TRUE, TRUE, Op.Max_Pos}
 ,{[4,5,6,7], TRUE, FALSE, Op.Max_Pos}
 ,{[0,4,8], TRUE, TRUE, Op.Min_Pos}
 ,{[8,9,10,11], TRUE, FALSE, Op.Min_Pos}
 ,{[10/4, 26/4, 1/4], TRUE, TRUE, Op.Ave_Val}
 ,{[6/3, 8/3, 10/3, 13/3], TRUE, FALSE, Op.Ave_Val}
 ,{[10, 26, 1], TRUE, TRUE, Op.Sum_Val}
 ,{[6, 8, 10, 13], TRUE, FALSE, Op.Sum_Val}
 ,{[10/4, 26/4, 1/4], FALSE, TRUE, Op.Ave_Abs}
 ,{[6/3, 8/3, 10/3, 13/3], FALSE, FALSE, Op.Ave_Abs}
 ,{[10, 26, 1], FALSE, TRUE, Op.Sum_Abs}
 ,{[6, 8, 10, 13], FALSE, FALSE, Op.Sum_Abs}
 ,{[30, 174, 1], TRUE, TRUE, Op.Sum_Sqs}
 ,{[26, 40, 58, 81], TRUE, FALSE, Op.Sum_Sqs}
 ,{[0, 0, 3], TRUE, TRUE, Op.Zeros}
 ,{[1, 1, 1, 0], TRUE, FALSE, Op.Zeros}
 ,{[4, 4, 1], TRUE, TRUE, Op.Non_Zeros}
 ,{[2, 2, 2, 3], TRUE, FALSE, Op.Non_Zeros}], Test_In);

t_run := PROJECT(test_set, ex(LEFT));

EXPORT unit_test_gemagg := t_run;