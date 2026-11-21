(* tests/proof/angle.v: the reductions angles held in degrees make without error, proved in Rocq
   with Flocq.

   On the binary64 grid (53 bits, subnormals down to 2^-1074, round to nearest even), each step
   below is the C++ of source/direction.cpp, source/heading.cpp and source/turning.hpp, operation
   for operation, and each theorem says the step is exact: what the code computes is the real
   value written, not a rounding of it.

     Quadrant  rest = high - 90 * rint(high * RN(1/90))       direction.cpp, |high| < 2^45
     Split     part = rest - rint(rest * 8) * 0.125           direction.cpp
     Arc       top  = t - rint(t * 128) * 2^-7                heading.cpp, 0 <= t <= 1
     Wrap      w - 360 for 180 < w <= 540, w + 360 for -540 < w <= -180   turning.hpp

   It also bounds what the steps leave: rest within 45 + 1/16, so the table index |k| <= 360, and
   part within 1/16 and t's step within 1/256, so the polynomials' remainders hold.

   One lemma does the work, Aligned: a double x of at least 2^a is a whole multiple of 2^(a-52);
   take from it any multiple of 2^b, b >= a - 52, and if the difference is below 2^(a+1) it has
   at most 53 bits, so it is a double. *)

From Stdlib Require Import Reals Lra Lia ZArith.
From Flocq Require Import Core Relative Mult_error.

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

(* std::rint under round to nearest even, as an integer. *)
Definition rint (x : R) : Z := ZnearestE x.

(* The powers of two the steps use. *)
Lemma two_1 : bpow radix2 1 = 2.     Proof. reflexivity. Qed.
Lemma two_3 : bpow radix2 3 = 8.     Proof. reflexivity. Qed.
Lemma two_5 : bpow radix2 5 = 32.    Proof. reflexivity. Qed.
Lemma two_6 : bpow radix2 6 = 64.    Proof. reflexivity. Qed.
Lemma two_7 : bpow radix2 7 = 128.   Proof. reflexivity. Qed.
Lemma two_8 : bpow radix2 8 = 256.   Proof. reflexivity. Qed.
Lemma two_m3 : bpow radix2 (-3) = / 8.   Proof. reflexivity. Qed.
Lemma two_m4 : bpow radix2 (-4) = / 16.  Proof. reflexivity. Qed.
Lemma two_m7 : bpow radix2 (-7) = / 128. Proof. reflexivity. Qed.
Lemma two_m8 : bpow radix2 (-8) = / 256. Proof. reflexivity. Qed.

(* A double of at least 2^a is a whole multiple of 2^(a-52). *)
Lemma Multiple :
  forall x a,
  fmt x -> (emin <= a - 52)%Z -> bpow radix2 a <= Rabs x ->
  exists m, x = IZR m * bpow radix2 (a - 52).
Proof.
  intros x a Fx Ha Hx.
  assert (Hm : (a + 1 <= mag radix2 x)%Z).
  { apply mag_ge_bpow. now replace (a + 1 - 1)%Z with a by ring. }
  assert (He : (a - 52 <= cexp radix2 fexp x)%Z).
  { unfold cexp, fexp, FLT_exp, prec. lia. }
  exists (Ztrunc (scaled_mantissa radix2 fexp x)
          * Zpower radix2 (cexp radix2 fexp x - (a - 52)))%Z.
  unfold fmt, generic_format in Fx.
  rewrite Fx at 1. unfold F2R; simpl.
  rewrite mult_IZR, (IZR_Zpower radix2) by lia.
  rewrite Rmult_assoc, <- bpow_plus.
  do 2 f_equal. ring.
Qed.

(* The difference of a double of at least 2^a and a multiple of 2^b, b >= a - 52, is a double
   when it is below 2^(a+1). *)
Lemma Aligned :
  forall x a b k,
  fmt x -> (emin <= a - 52)%Z -> (a - 52 <= b)%Z ->
  bpow radix2 a <= Rabs x ->
  Rabs (x - IZR k * bpow radix2 b) < bpow radix2 (a + 1) ->
  fmt (x - IZR k * bpow radix2 b).
Proof.
  intros x a b k Fx Ha Hb Hx Hd.
  destruct (Multiple x a Fx Ha Hx) as [m Hm].
  set (j := (k * Zpower radix2 (b - (a - 52)))%Z).
  assert (Hy : IZR k * bpow radix2 b = IZR j * bpow radix2 (a - 52)).
  { unfold j. rewrite mult_IZR, (IZR_Zpower radix2) by lia.
    rewrite Rmult_assoc, <- bpow_plus. do 2 f_equal. ring. }
  assert (Hxy : x - IZR k * bpow radix2 b = IZR (m - j) * bpow radix2 (a - 52)).
  { rewrite Hy, minus_IZR. rewrite Hm at 1. ring. }
  rewrite Hxy in Hd |- *.
  assert (Hl : Rabs (IZR (m - j)) < bpow radix2 prec).
  { rewrite Rabs_mult, (Rabs_pos_eq (bpow radix2 (a - 52))) in Hd by apply bpow_ge_0.
    replace (a + 1)%Z with (prec + (a - 52))%Z in Hd by (unfold prec; ring).
    rewrite bpow_plus in Hd.
    apply Rmult_lt_reg_r with (bpow radix2 (a - 52)); [apply bpow_gt_0 | exact Hd]. }
  rewrite <- abs_IZR, <- (IZR_Zpower radix2) in Hl by (unfold prec; lia).
  apply lt_IZR in Hl.
  apply generic_format_FLT.
  apply FLT_spec with (Float radix2 (m - j) (a - 52)).
  - reflexivity.
  - exact Hl.
  - exact Ha.
Qed.

(* A whole number below 2^53 times a power of two is a double. *)
Lemma Whole :
  forall k e, (Z.abs k < 2 ^ 53)%Z -> (emin <= e)%Z -> fmt (IZR k * bpow radix2 e).
Proof.
  intros k e Hk He.
  apply generic_format_FLT.
  apply FLT_spec with (Float radix2 k e); simpl; [reflexivity | | exact He].
  unfold prec. exact Hk.
Qed.

Lemma Exact : forall x, fmt x -> rnd x = x.
Proof. intros x Fx. now apply round_generic; auto with typeclass_instances. Qed.

(* rint picks the nearest whole number. *)
Lemma Nearest : forall x, Rabs (x - IZR (rint x)) <= / 2.
Proof. intros x. apply Znearest_half. Qed.

Lemma Nearest_zero : forall x, Rabs x < / 2 -> rint x = 0%Z.
Proof.
  intros x Hx. apply Znearest_imp. now rewrite Rminus_0_r.
Qed.

(* ---------------------------------------------------------------------------------------------
   Wrap: turning.hpp's Wrap steps back by a turn exactly. *)

Theorem Wrap_high : forall w, fmt w -> 180 < w <= 540 -> fmt (w - 360).
Proof.
  intros w Fw Hw.
  replace (w - 360) with (w - IZR 45 * bpow radix2 3) by (rewrite two_3; ring).
  assert (G1 : bpow radix2 7 <= Rabs w) by (rewrite two_7, Rabs_pos_eq; lra).
  assert (G2 : Rabs (w - IZR 45 * bpow radix2 3) < bpow radix2 (7 + 1)).
  { replace (7 + 1)%Z with 8%Z by ring. rewrite two_8, two_3. apply Rabs_def1; lra. }
  exact (Aligned w 7 3 45 Fw ltac:(unfold emin; lia) ltac:(lia) G1 G2).
Qed.

Theorem Wrap_low : forall w, fmt w -> -540 < w <= -180 -> fmt (w + 360).
Proof.
  intros w Fw Hw.
  replace (w + 360) with (w - IZR (-45) * bpow radix2 3) by (rewrite two_3; ring).
  assert (G1 : bpow radix2 7 <= Rabs w) by (rewrite two_7, Rabs_left; lra).
  assert (G2 : Rabs (w - IZR (-45) * bpow radix2 3) < bpow radix2 (7 + 1)).
  { replace (7 + 1)%Z with 8%Z by ring. rewrite two_8, two_3. apply Rabs_def1; lra. }
  exact (Aligned w 7 3 (-45) Fw ltac:(unfold emin; lia) ltac:(lia) G1 G2).
Qed.

Print Assumptions Wrap_high.
Print Assumptions Wrap_low.

(* ---------------------------------------------------------------------------------------------
   Split: direction.cpp's step from the quarter-turn rest to an eighth of a degree. *)

Theorem Split :
  forall r, fmt r -> Rabs r < 45 + / 16 ->
  let k := rint (rnd (r * 8)) in
  rnd (r - rnd (IZR k * / 8)) = r - IZR k / 8
  /\ Rabs (r - IZR k / 8) <= / 16
  /\ (Z.abs k <= 360)%Z.
Proof.
  intros r Fr Hr k.
  assert (H8 : rnd (r * 8) = r * 8).
  { apply Exact. rewrite <- two_3. apply mult_bpow_pos_exact_FLT; auto. lia. }
  assert (Hk : Rabs (r * 8 - IZR k) <= / 2).
  { unfold k. rewrite H8. apply Nearest. }
  assert (Hd : Rabs (r - IZR k / 8) <= / 16).
  { apply Rabs_le_inv in Hk. apply Rabs_le. lra. }
  assert (Hb : (Z.abs k <= 360)%Z).
  { apply Rabs_le_inv in Hk. apply Rabs_lt_inv in Hr.
    assert (Ha : (-361 < k < 361)%Z) by (split; apply lt_IZR; lra).
    lia. }
  assert (Hq : rnd (IZR k * / 8) = IZR k / 8).
  { unfold Rdiv. apply Exact. rewrite <- two_m3. apply Whole; [lia | unfold emin; lia]. }
  split; [|split; assumption].
  rewrite Hq. apply Exact.
  destruct (Rlt_or_le (Rabs r) (/ 16)) as [Hs | Hl].
  - (* below 1/16: rint (8 r) is 0, and the part is r itself *)
    assert (Hz : k = 0%Z).
    { unfold k. rewrite H8. apply Nearest_zero.
      apply Rabs_lt_inv in Hs. apply Rabs_lt. lra. }
    rewrite Hz. replace (r - IZR 0 / 8) with r by (simpl; field). exact Fr.
  - replace (r - IZR k / 8) with (r - IZR k * bpow radix2 (-3)) by (rewrite two_m3; field).
    assert (G1 : bpow radix2 (-4) <= Rabs r) by (rewrite two_m4; lra).
    assert (G2 : Rabs (r - IZR k * bpow radix2 (-3)) < bpow radix2 (-4 + 1)).
    { replace (-4 + 1)%Z with (-3)%Z by ring. rewrite two_m3.
      replace (r - IZR k * / 8) with (r - IZR k / 8) by field. lra. }
    exact (Aligned r (-4) (-3) k Fr ltac:(unfold emin; lia) ltac:(lia) G1 G2).
Qed.

Print Assumptions Split.

(* ---------------------------------------------------------------------------------------------
   Arc: heading.cpp's step from the ratio to the nearest 1/128. *)

Theorem Arc :
  forall t, fmt t -> 0 <= t <= 1 ->
  let k := rint (rnd (t * 128)) in
  rnd (t - rnd (IZR k * / 128)) = t - IZR k / 128
  /\ Rabs (t - IZR k / 128) <= / 256
  /\ (0 <= k <= 128)%Z.
Proof.
  intros t Ft Ht k.
  assert (H8 : rnd (t * 128) = t * 128).
  { apply Exact. rewrite <- two_7. apply mult_bpow_pos_exact_FLT; auto. lia. }
  assert (Hk : Rabs (t * 128 - IZR k) <= / 2).
  { unfold k. rewrite H8. apply Nearest. }
  assert (Hd : Rabs (t - IZR k / 128) <= / 256).
  { apply Rabs_le_inv in Hk. apply Rabs_le. lra. }
  assert (Hb : (0 <= k <= 128)%Z).
  { apply Rabs_le_inv in Hk.
    assert (Ha : (-1 < k < 129)%Z) by (split; apply lt_IZR; simpl; lra).
    lia. }
  assert (Hq : rnd (IZR k * / 128) = IZR k / 128).
  { unfold Rdiv. apply Exact. rewrite <- two_m7. apply Whole; [lia | unfold emin; lia]. }
  split; [|split; assumption].
  rewrite Hq. apply Exact.
  destruct (Rlt_or_le t (/ 256)) as [Hs | Hl].
  - assert (Hz : k = 0%Z).
    { unfold k. rewrite H8. apply Nearest_zero. apply Rabs_lt. lra. }
    rewrite Hz. replace (t - IZR 0 / 128) with t by (simpl; field). exact Ft.
  - replace (t - IZR k / 128) with (t - IZR k * bpow radix2 (-7)) by (rewrite two_m7; field).
    assert (G1 : bpow radix2 (-8) <= Rabs t) by (rewrite two_m8, Rabs_pos_eq; lra).
    assert (G2 : Rabs (t - IZR k * bpow radix2 (-7)) < bpow radix2 (-8 + 1)).
    { replace (-8 + 1)%Z with (-7)%Z by ring. rewrite two_m7.
      replace (t - IZR k * / 128) with (t - IZR k / 128) by field. lra. }
    exact (Aligned t (-8) (-7) k Ft ltac:(unfold emin; lia) ltac:(lia) G1 G2).
Qed.

Print Assumptions Arc.

(* ---------------------------------------------------------------------------------------------
   Quadrant: direction.cpp's step from any angle below 2^45 degrees to within 45 + 1/16 of 0. The
   quarter turn is picked with a rounded 1/90 and a rounded product, so this one bounds rounding
   errors before the alignment argument applies. *)

Lemma Ninetieth :
  let c := rnd (/ 90) in 0 < c /\ Rabs (c - / 90) <= / 2 * bpow radix2 (-52) * / 90.
Proof.
  intros c.
  assert (Hn : bpow radix2 (emin + prec - 1) <= Rabs (/ 90)).
  { rewrite Rabs_pos_eq by lra.
    apply Rle_trans with (bpow radix2 (-7)).
    - apply bpow_le. unfold emin, prec. lia.
    - rewrite two_m7. lra. }
  pose proof (relative_error_N_FLT radix2 emin prec prec_gt_0 (fun x => negb (Z.even x)) (/ 90) Hn)
    as He.
  fold fexp in He. fold rnd in He. fold c in He.
  replace (- prec + 1)%Z with (-52)%Z in He by (unfold prec; ring).
  rewrite (Rabs_pos_eq (/ 90)) in He by lra.
  split; [|exact He].
  assert (Hs : bpow radix2 (-52) <= 1).
  { replace 1 with (bpow radix2 0) by reflexivity. apply bpow_le. lia. }
  apply Rabs_le_inv in He. lra.
Qed.

Theorem Quadrant :
  forall h, fmt h -> Rabs h < bpow radix2 45 ->
  let q := rint (rnd (h * rnd (/ 90))) in
  rnd (h - rnd (90 * IZR q)) = h - 90 * IZR q
  /\ Rabs (h - 90 * IZR q) < 45 + / 16.
Proof.
  intros h Fh Hh q.
  destruct Ninetieth as [Hc0 Hc].
  set (c := rnd (/ 90)) in *.
  set (p := rnd (h * c)).
  assert (Hu : bpow radix2 (-52) = / 4503599627370496) by reflexivity.
  assert (Ht : bpow radix2 45 = 35184372088832) by reflexivity.
  rewrite Hu in Hc. rewrite Ht in Hh.
  (* p = h c (1 + eps) + eta *)
  destruct (error_N_FLT radix2 emin prec prec_gt_0 (fun x => negb (Z.even x)) (h * c))
    as [eps [eta [He [Hn [_ Hp]]]]].
  fold fexp rnd p in Hp.
  replace (- prec + 1)%Z with (-52)%Z in He by (unfold prec; ring).
  rewrite Hu in He.
  assert (Hm : bpow radix2 emin <= / 1152921504606846976).
  { replace (/ 1152921504606846976) with (bpow radix2 (-60)) by reflexivity.
    apply bpow_le. unfold emin. lia. }
  assert (Hm0 : 0 < bpow radix2 emin) by apply bpow_gt_0.
  set (tau := bpow radix2 emin) in *.
  assert (Hn' : Rabs eta <= tau) by lra.
  (* |h| |c - 1/90| and |h c| bounded *)
  assert (Ha : Rabs h * Rabs (c - / 90) <= 35184372088832 * (/ 2 * / 4503599627370496 * / 90)).
  { apply Rmult_le_compat; try apply Rabs_pos; lra. }
  assert (Hcb : c <= / 90 + / 2 * / 4503599627370496 * / 90).
  { apply Rabs_le_inv in Hc. lra. }
  assert (Hhc : Rabs (h * c) <= Rabs h * (/ 90 + / 2 * / 4503599627370496 * / 90)).
  { rewrite Rabs_mult, (Rabs_pos_eq c) by lra. apply Rmult_le_compat_l; [apply Rabs_pos | lra]. }
  assert (Hpe : Rabs (p - h * c) <= Rabs (h * c) * / 2 * / 4503599627370496 + tau).
  { rewrite Hp. replace (h * c * (1 + eps) + eta - h * c) with (h * c * eps + eta) by ring.
    eapply Rle_trans; [apply Rabs_triang|]. rewrite Rabs_mult.
    apply Rplus_le_compat; [|exact Hn'].
    rewrite Rmult_assoc. apply Rmult_le_compat_l; [apply Rabs_pos | lra]. }
  assert (Hq : Rabs (p - IZR q) <= / 2) by (unfold q; apply Nearest).
  (* h - 90 q = 90 (h (1/90 - c) + (h c - p) + (p - q)) *)
  assert (Hb : Rabs (h - 90 * IZR q) < 45 + / 16).
  { replace (h - 90 * IZR q) with (90 * (h * (/ 90 - c) + (h * c - p) + (p - IZR q))) by field.
    rewrite Rabs_mult, (Rabs_pos_eq 90) by lra.
    assert (H3 : Rabs (h * (/ 90 - c) + (h * c - p) + (p - IZR q))
                 <= Rabs h * Rabs (c - / 90) + Rabs (p - h * c) + Rabs (p - IZR q)).
    { eapply Rle_trans; [apply Rabs_triang|]. apply Rplus_le_compat_r.
      eapply Rle_trans; [apply Rabs_triang|].
      rewrite Rabs_mult, (Rabs_minus_sym (/ 90)), (Rabs_minus_sym (h * c)). lra. }
    lra. }
  split; [|exact Hb].
  (* 90 q is a double: |q| is below 2^45 *)
  assert (Hq90 : rnd (90 * IZR q) = 90 * IZR q).
  { apply Exact.
    replace (90 * IZR q) with (IZR (45 * q) * bpow radix2 1) by (rewrite two_1, mult_IZR; ring).
    apply Whole; [|unfold emin; lia].
    assert (Hz : Rabs (IZR (45 * q)) < 2 ^ 53).
    { rewrite mult_IZR, Rabs_mult, (Rabs_pos_eq (IZR 45)) by (apply IZR_le; lia).
      assert (Hr : Rabs (90 * IZR q) < 35184372088832 + 46).
      { replace (90 * IZR q) with (h - (h - 90 * IZR q)) by ring.
        eapply Rle_lt_trans; [apply Rabs_triang|]. rewrite Rabs_Ropp. lra. }
      rewrite Rabs_mult, (Rabs_pos_eq 90) in Hr by lra. simpl. lra. }
    rewrite <- abs_IZR in Hz. replace (2 ^ 53) with (IZR (2 ^ 53)) in Hz by (simpl; lra).
    apply lt_IZR in Hz. exact Hz. }
  rewrite Hq90. apply Exact.
  destruct (Rlt_or_le (Rabs h) 32) as [Hs | Hl].
  - (* below 32 degrees: the quarter is 0 *)
    assert (Hz : q = 0%Z).
    { unfold q. apply Nearest_zero. fold c p.
      replace p with (h * c + (p - h * c)) by ring.
      eapply Rle_lt_trans; [apply Rabs_triang|].
      assert (Hs' : Rabs (h * c) <= 32 * (/ 90 + / 2 * / 4503599627370496 * / 90)).
      { eapply Rle_trans; [exact Hhc|]. apply Rmult_le_compat_r; lra. }
      lra. }
    rewrite Hz. replace (h - 90 * IZR 0) with h by (simpl; ring). exact Fh.
  - replace (h - 90 * IZR q) with (h - IZR (45 * q) * bpow radix2 1)
      by (rewrite two_1, mult_IZR; ring).
    assert (G1 : bpow radix2 5 <= Rabs h) by (rewrite two_5; lra).
    assert (G2 : Rabs (h - IZR (45 * q) * bpow radix2 1) < bpow radix2 (5 + 1)).
    { replace (5 + 1)%Z with 6%Z by ring. rewrite two_6, two_1.
      replace (h - IZR (45 * q) * 2) with (h - 90 * IZR q) by (rewrite mult_IZR; ring). lra. }
    exact (Aligned h 5 1 (45 * q) Fh ltac:(unfold emin; lia) ltac:(lia) G1 G2).
Qed.

Print Assumptions Quadrant.
