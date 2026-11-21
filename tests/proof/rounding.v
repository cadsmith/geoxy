(* tests/proof/rounding.v: the IEEE 754 binary64 facts the exact predicates' proofs rest on,
   proved in Rocq with Flocq. tests/proof/lemma.py states them to Z3 as axioms; here each is a
   theorem about round to nearest, ties to even, on the binary64 grid: 53 bits of precision,
   subnormals down to 2^-1074, exponents unbounded above, overflow being left to Z3's IEEE 754
   theory. The error-free transformations are source/exact/two.hpp's, operation for operation:
   Flocq proves Knuth's TwoSum and Dekker's product with Veltkamp's split in its own order, and each theorem here rewrites the code's order into Flocq's. *)

From Stdlib Require Import Reals Lra Lia ZArith.
From Flocq Require Import Core Plus_error Mult_error Relative Pff2Flocq.
From Flocq Require Import IEEE754.BinarySingleNaN IEEE754.Binary.

Open Scope R_scope.

Definition emin : Z := (-1074)%Z.
Definition prec : Z := 53%Z.

#[global] Instance prec_gt_0 : Prec_gt_0 prec.
Proof. now unfold Prec_gt_0, prec. Qed.

Definition fexp := FLT_exp emin prec.

#[global] Instance fexp_valid : Valid_exp fexp.
Proof. unfold fexp; apply FLT_exp_valid; auto with typeclass_instances. Qed.
Definition fmt := generic_format radix2 fexp.
Definition rnd := round radix2 fexp ZnearestE.

(* u, half the distance from 1 to the next double; uN, the smallest normal double; h, half the
   smallest subnormal. *)
Definition u := bpow radix2 (-53).
Definition uN := bpow radix2 (-1022).
Definition h := bpow radix2 (-1075).

Lemma half_bpow : forall e, / 2 * bpow radix2 e = bpow radix2 (e - 1).
Proof.
  intros e; replace (e - 1)%Z with (e + - (1))%Z by lia.
  rewrite bpow_plus, bpow_opp, bpow_1; simpl (IZR (radix_val radix2)); lra.
Qed.

Lemma half_ulp_one : / 2 * bpow radix2 (- prec + 1) = u.
Proof. rewrite half_bpow; unfold u, prec; f_equal. Qed.

Lemma choice_sym :
  forall x, negb (Z.even x) = negb (negb (Z.even (- (x + 1)))).
Proof.
  intros x; rewrite Z.even_opp, Z.even_add; simpl.
  destruct (Z.even x); reflexivity.
Qed.

Lemma format_uN : fmt uN.
Proof.
  apply generic_format_bpow; unfold fexp, FLT_exp, emin, prec; lia.
Qed.

Lemma round_zero : rnd 0 = 0.
Proof. apply round_0; auto with typeclass_instances. Qed.

(* The relative error of a rounding clear of underflow, against the exact value and the result:
   lemma.py's rounded, for products. *)
Theorem rounded_relative :
  forall z, (z = 0 \/ uN <= Rabs z) ->
  Rabs (rnd z - z) <= u * Rabs z /\ Rabs (rnd z - z) <= u * Rabs (rnd z).
Proof.
  intros z [Hz | Hz].
  - subst z; rewrite round_zero, Rminus_0_r, Rabs_R0; lra.
  - assert (H : bpow radix2 (emin + prec - 1) <= Rabs z) by exact Hz.
    split.
    + rewrite <- half_ulp_one; now apply relative_error_N_FLT.
    + rewrite <- half_ulp_one; now apply relative_error_N_FLT_round.
Qed.

(* A sum or difference of two doubles: where it would underflow it is exact, so both bounds hold
   always. lemma.py's rounded, for sums. *)
Theorem sum_relative :
  forall x y, fmt x -> fmt y ->
  Rabs (rnd (x + y) - (x + y)) <= u * Rabs (x + y) /\
  Rabs (rnd (x + y) - (x + y)) <= u * Rabs (rnd (x + y)).
Proof.
  intros x y Fx Fy.
  destruct (Rlt_or_le (Rabs (x + y)) uN) as [Hs | Hs].
  - assert (F : fmt (x + y)).
    { apply (FLT_format_plus_small radix2 emin prec); auto.
      apply Rlt_le, Rlt_le_trans with (1 := Hs).
      apply bpow_le; unfold emin, prec; lia. }
    unfold rnd; rewrite round_generic; auto with typeclass_instances.
    rewrite Rminus_diag, Rabs_R0; split; apply Rmult_le_pos; try apply Rabs_pos;
      apply bpow_ge_0.
  - apply rounded_relative; now right.
Qed.

(* A rounding of any real, gradual underflow included: lemma.py's multiplied. *)
Theorem multiplied :
  forall z,
  (uN <= Rabs z ->
   Rabs (rnd z - z) <= u * Rabs z /\ Rabs (rnd z - z) <= u * Rabs (rnd z)) /\
  (Rabs z < uN -> Rabs (rnd z - z) <= h /\ Rabs (rnd z) <= uN) /\
  (0 <= z -> 0 <= rnd z) /\
  (z <= 0 -> rnd z <= 0).
Proof.
  intros z; split; [| split; [| split]].
  - intros Hz; apply rounded_relative; now right.
  - intros Hz; split.
    + apply Rle_trans with (1 := error_le_half_ulp radix2 fexp _ z).
      unfold fexp; rewrite ulp_FLT_small; auto with typeclass_instances.
      * right; rewrite half_bpow; unfold h, emin; f_equal.
      * apply Rlt_le_trans with (1 := Hz); apply bpow_le; unfold emin, prec; lia.
    + destruct (Rabs_def2 _ _ Hz) as [Hlo Hhi]; apply Rabs_le; split.
      * rewrite <- (Ropp_involutive (rnd z)); apply Ropp_le_contravar.
        unfold rnd; rewrite <- round_NE_opp.
        apply round_le_generic; auto with typeclass_instances; [apply format_uN | lra].
      * unfold rnd.
        apply round_le_generic; auto with typeclass_instances; [apply format_uN | lra].
  - intros Hz; unfold rnd.
    apply round_ge_generic; auto with typeclass_instances; apply generic_format_0.
  - intros Hz; unfold rnd.
    apply round_le_generic; auto with typeclass_instances; apply generic_format_0.
Qed.

(* Rounding is monotone and odd. *)
Theorem monotone : forall z1 z2, z1 <= z2 -> rnd z1 <= rnd z2.
Proof. intros; apply round_le; auto with typeclass_instances. Qed.

Theorem symmetric : forall z, rnd (- z) = - rnd z.
Proof. intros; apply round_NE_opp. Qed.

Lemma rnd_ext : forall a b, a = b -> rnd a = rnd b.
Proof. intros; now subst. Qed.

(* TwoSum, two.hpp's order: sum, virtualSecond, virtualFirst, secondRound, firstRound. *)
Theorem two_sum :
  forall a b, fmt a -> fmt b ->
  let s := rnd (a + b) in
  let vb := rnd (s - a) in
  let va := rnd (s - vb) in
  let rb := rnd (b - vb) in
  let ra := rnd (a - va) in
  s + rnd (ra + rb) = a + b.
Proof.
  intros a b Fa Fb; simpl.
  apply (TwoSum_correct emin prec (fun x => negb (Z.even x))); auto.
  - unfold prec; lia.
  - unfold emin; lia.
  - apply choice_sym.
Qed.

(* TwoDiffTail, two.hpp's order: virtualSecond, virtualFirst, secondRound, firstRound; TwoSum
   of a and -b, its operations negated where rounding is odd. *)
Theorem two_diff_tail :
  forall a b, fmt a -> fmt b ->
  let d := rnd (a - b) in
  let vs := rnd (a - d) in
  let vf := rnd (d + vs) in
  let sr := rnd (vs - b) in
  let fr := rnd (a - vf) in
  d + rnd (fr + sr) = a - b.
Proof.
  intros a b Fa Fb; simpl.
  assert (Fb' : fmt (- b)) by now apply generic_format_opp.
  pose proof (TwoSum_correct emin prec (fun x => negb (Z.even x)) ltac:(unfold prec; lia)
                ltac:(unfold emin; lia) choice_sym a (- b) Fa Fb') as T.
  change (rnd (a + - b) + rnd (rnd (a - rnd (rnd (a + - b) - rnd (rnd (a + - b) - a)))
          + rnd (- b - rnd (rnd (a + - b) - a))) = a + - b) in T.
  replace (a + - b) with (a - b) in T by ring.
  assert (Hx : rnd (rnd (a - b) - a) = - rnd (a - rnd (a - b))).
  { rewrite <- symmetric; apply rnd_ext; ring. }
  rewrite Hx in T.
  replace (rnd (a - b) - - rnd (a - rnd (a - b)))
    with (rnd (a - b) + rnd (a - rnd (a - b))) in T by ring.
  replace (- b - - rnd (a - rnd (a - b)))
    with (rnd (a - rnd (a - b)) - b) in T by ring.
  exact T.
Qed.

(* FastTwoSum, twofold.hpp's order: the sum, then small - (sum - big), exact where the small term
   is no larger than the big one (Dekker). Flocq's takes small + rnd (big - sum), the same value,
   rounding being odd. *)
Theorem fast_two_sum :
  forall a b, fmt a -> fmt b -> Rabs b <= Rabs a ->
  let s := rnd (a + b) in
  s + rnd (b - rnd (s - a)) = a + b.
Proof.
  intros a b Fa Fb Hab; simpl.
  pose proof (Fast2Sum_correct emin prec (fun x => negb (Z.even x)) ltac:(unfold prec; lia)
                ltac:(unfold emin; lia) choice_sym a b Fa Fb Hab) as T.
  assert (Hx : rnd (rnd (a + b) - a) = - rnd (a - rnd (a + b))).
  { rewrite <- symmetric; apply rnd_ext; ring. }
  rewrite Hx.
  replace (b - - rnd (a - rnd (a + b))) with (b + rnd (a - rnd (a + b))) by ring.
  exact T.
Qed.

(* TwoProduct where the hardware fuses: fma(a, b, -p) rounds a b - p, a double. *)
Theorem two_product_fma :
  forall a b, fmt a -> fmt b ->
  (a * b <> 0 -> bpow radix2 (-969) <= Rabs (a * b)) ->
  let p := rnd (a * b) in
  p + rnd (a * b - p) = a * b.
Proof.
  intros a b Fa Fb Hab; simpl.
  assert (F : fmt (rnd (a * b) - a * b)).
  { apply mult_error_FLT; auto with typeclass_instances. }
  assert (F' : fmt (a * b - rnd (a * b))).
  { replace (a * b - rnd (a * b)) with (- (rnd (a * b) - a * b)) by ring.
    now apply generic_format_opp. }
  unfold rnd at 2; rewrite round_generic; auto with typeclass_instances; ring.
Qed.

(* orient.cpp's Safe: a rounded product at least 2^-968 in magnitude is one at least 2^-969 before
   rounding, which two_product_fma asks. *)
Theorem safe_product :
  forall z, bpow radix2 (-968) <= Rabs (rnd z) -> bpow radix2 (-969) <= Rabs z.
Proof.
  intros z Hz.
  destruct (Rle_or_lt (bpow radix2 (-969)) (Rabs z)) as [H | H]; [exact H | exfalso].
  assert (Hr : Rabs (rnd z) <= bpow radix2 (-969)).
  { unfold rnd; rewrite <- round_NE_abs; auto with typeclass_instances.
    apply round_le_generic; auto with typeclass_instances; [|lra].
    apply generic_format_bpow; unfold fexp, FLT_exp, emin, prec; lia. }
  assert (bpow radix2 (-969) < bpow radix2 (-968)) by (apply bpow_lt; lia).
  lra.
Qed.

(* orient.cpp's CrossSign: with p = c d rounded, e = fma(c, d, -p) and f = fma(a, b, -p), f > e
   proves a b - c d positive and f < e negative. e is c d - p exactly; were a b - p at most e,
   f, its rounding, would be at most e, by monotony. *)
Theorem cross_sign :
  forall a b c d, fmt a -> fmt b -> fmt c -> fmt d ->
  (c * d <> 0 -> bpow radix2 (-969) <= Rabs (c * d)) ->
  let p := rnd (c * d) in
  let e := rnd (c * d - p) in
  let f := rnd (a * b - p) in
  (f > e -> a * b - c * d > 0) /\ (f < e -> a * b - c * d < 0).
Proof.
  intros a b c d Fa Fb Fc Fd Hcd p e f.
  assert (Fe : fmt (c * d - p)).
  { unfold p; replace (c * d - rnd (c * d)) with (- (rnd (c * d) - c * d)) by ring.
    apply generic_format_opp; apply mult_error_FLT; auto with typeclass_instances. }
  assert (He : e = c * d - p).
  { unfold e, rnd; apply round_generic; auto with typeclass_instances. }
  split; intros Hf.
  - destruct (Rle_or_lt (a * b - p) e) as [H | H]; [|lra].
    exfalso; assert (f <= rnd e) by (apply monotone; exact H).
    assert (rnd e = e) by (unfold rnd; rewrite He; apply round_generic;
                           auto with typeclass_instances).
    lra.
  - destruct (Rle_or_lt e (a * b - p)) as [H | H]; [|lra].
    exfalso; assert (rnd e <= f) by (apply monotone; exact H).
    assert (rnd e = e) by (unfold rnd; rewrite He; apply round_generic;
                           auto with typeclass_instances).
    lra.
Qed.

(* Veltkamp's split, as Flocq orders it, and as two.hpp does: scaled, big, high, low, where
   Flocq's q is -big. *)
Definition splitter := bpow radix2 27 + 1.
Definition high v := rnd (rnd (v - rnd (v * splitter)) + rnd (v * splitter)).
Definition low v := rnd (v - high v).

Theorem split_code :
  forall v,
  let scaled := rnd (splitter * v) in
  let big := rnd (scaled - v) in
  high v = rnd (scaled - big) /\ low v = rnd (v - rnd (scaled - big)).
Proof.
  intros v; simpl.
  assert (H : high v = rnd (rnd (splitter * v) - rnd (rnd (splitter * v) - v))).
  { unfold high; replace (v * splitter) with (splitter * v) by ring.
    apply rnd_ext.
    replace (v - rnd (splitter * v)) with (- (rnd (splitter * v) - v)) by ring.
    rewrite symmetric; ring. }
  split; [exact H | unfold low; now rewrite H].
Qed.

(* TwoProduct where products are split, two.hpp's order: error1, error2, error3, small. It is
   Flocq's Dekker with the factors swapped and its partial sums negated. *)
Theorem two_product_split :
  forall a b, fmt a -> fmt b ->
  (a * b = 0 \/ bpow radix2 (-969) <= Rabs (a * b)) ->
  let p := rnd (a * b) in
  let e1 := rnd (p - rnd (high a * high b)) in
  let e2 := rnd (e1 - rnd (low a * high b)) in
  let e3 := rnd (e2 - rnd (high a * low b)) in
  p + rnd (rnd (low a * low b) - e3) = a * b.
Proof.
  intros a b Fa Fb Hab; simpl.
  rewrite (Rmult_comm a b) in Hab.
  pose proof (Dekker radix2 emin prec (fun x => negb (Z.even x))
                ltac:(unfold prec; lia) ltac:(unfold emin; lia) b a Fb Fa
                (or_introl eq_refl)) as [D _].
  specialize (D Hab).
  change (b * a = rnd (b * a) + rnd (rnd (rnd (rnd (- rnd (b * a) + rnd (high b * high a))
          + rnd (high b * low a)) + rnd (low b * high a)) + rnd (low b * low a))) in D.
  rewrite (Rmult_comm b a), (Rmult_comm (high b) (high a)), (Rmult_comm (high b) (low a)),
          (Rmult_comm (low b) (high a)), (Rmult_comm (low b) (low a)) in D.
  set (r := rnd (a * b)) in *.
  set (t1 := rnd (- r + rnd (high a * high b))) in *.
  set (t2 := rnd (t1 + rnd (low a * high b))) in *.
  set (t3 := rnd (t2 + rnd (high a * low b))) in *.
  assert (N : forall x y, rnd (- x - y) = - rnd (x + y)).
  { intros; rewrite <- symmetric; apply rnd_ext; ring. }
  assert (E1 : rnd (r - rnd (high a * high b)) = - t1).
  { unfold t1; rewrite <- symmetric; apply rnd_ext; ring. }
  rewrite E1.
  assert (E2 : rnd (- t1 - rnd (low a * high b)) = - t2) by (unfold t2; apply N).
  rewrite E2.
  assert (E3 : rnd (- t2 - rnd (high a * low b)) = - t3) by (unfold t3; apply N).
  rewrite E3.
  rewrite D.
  f_equal; apply rnd_ext; ring.
Qed.

(* IEEE 754 binary64 itself, Flocq's bit level model of it, signs, infinities and NaNs: an
   addition, a subtraction, a multiplication or a fused multiply-add of finite doubles that does
   not overflow is rnd of the exact value. So the theorems above hold of the hardware's
   operations. *)
Lemma prec_gt_53 : Prec_gt_0 53.
Proof. now unfold Prec_gt_0. Qed.

Lemma prec_lt_1024 : Prec_lt_emax 53 1024.
Proof. now unfold Prec_lt_emax. Qed.

Notation double := (binary_float 53 1024).
Notation value := (B2R 53 1024).
Notation finite := (is_finite 53 1024).

Theorem ieee_plus :
  forall nan (x y : double), finite x = true -> finite y = true ->
  Rabs (rnd (value x + value y)) < bpow radix2 1024 ->
  value (Bplus 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y) = rnd (value x + value y).
Proof.
  intros nan x y Fx Fy H.
  pose proof (Bplus_correct 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y Fx Fy) as C.
  rewrite Rlt_bool_true in C by exact H.
  now destruct C as [C _].
Qed.

Theorem ieee_minus :
  forall nan (x y : double), finite x = true -> finite y = true ->
  Rabs (rnd (value x - value y)) < bpow radix2 1024 ->
  value (Bminus 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y) = rnd (value x - value y).
Proof.
  intros nan x y Fx Fy H.
  pose proof (Bminus_correct 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y Fx Fy) as C.
  rewrite Rlt_bool_true in C by exact H.
  now destruct C as [C _].
Qed.

Theorem ieee_mult :
  forall nan (x y : double),
  Rabs (rnd (value x * value y)) < bpow radix2 1024 ->
  value (Bmult 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y) = rnd (value x * value y).
Proof.
  intros nan x y H.
  pose proof (Bmult_correct 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y) as C.
  rewrite Rlt_bool_true in C by exact H.
  now destruct C as [C _].
Qed.

Theorem ieee_fma :
  forall nan (x y z : double), finite x = true -> finite y = true -> finite z = true ->
  Rabs (rnd (value x * value y + value z)) < bpow radix2 1024 ->
  value (Bfma 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y z) =
    rnd (value x * value y + value z).
Proof.
  intros nan x y z Fx Fy Fz H.
  pose proof (Bfma_correct 53 1024 prec_gt_53 prec_lt_1024 nan mode_NE x y z Fx Fy Fz) as C.
  simpl in C.
  rewrite Rlt_bool_true in C by exact H.
  now destruct C as [C _].
Qed.

Theorem ieee_format : forall x : double, fmt (value x).
Proof. intros x; apply generic_format_B2R. Qed.

Print Assumptions rounded_relative.
Print Assumptions sum_relative.
Print Assumptions multiplied.
Print Assumptions monotone.
Print Assumptions two_sum.
Print Assumptions two_diff_tail.
Print Assumptions fast_two_sum.
Print Assumptions two_product_fma.
Print Assumptions safe_product.
Print Assumptions cross_sign.
Print Assumptions split_code.
Print Assumptions two_product_split.
Print Assumptions ieee_plus.
Print Assumptions ieee_minus.
Print Assumptions ieee_mult.
Print Assumptions ieee_fma.
Print Assumptions ieee_format.
