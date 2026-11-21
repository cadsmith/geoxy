(* tests/proof/integers.v: the whole-number arithmetic of source/exact/grid.hpp, grid.cpp and
   whole.cpp, proved in Rocq. Each definition here is one of theirs, step for step, in C++'s own
   semantics: std::uint64_t arithmetic is modulo 2^64, x << k is x 2^k modulo 2^64 and x >> k is
   x / 2^k for k below 64, a comparison is 1 or 0, ~x is 2^64 - 1 - x. A loop over words is a
   recursion over the list of them, lowest first, as the loop walks them. Each theorem says what
   the code computes is the true value; the last ones, that each of the four exact signs is the
   sign of the determinant of the doubles given, read from their bits as IEEE 754 binary64
   (Flocq's b64_of_bits). tests/proof/rocq.py builds this; tests/proof/integers.py runs the C++
   itself against exact rational arithmetic, many inputs, to tie the model to the code. *)

From Stdlib Require Import ZArith Lia List Bool Reals.
From Flocq Require Import Core IEEE754.Binary IEEE754.Bits.
Import ListNotations.

Open Scope Z_scope.

(* ---------------------------------------------------------------------------------------------- *)
(* Words.                                                                                         *)

Definition W : Z := 18446744073709551616.
Lemma W_pow : W = 2 ^ 64.
Proof. reflexivity. Qed.

Definition word x := 0 <= x < W.

Definition add x y := (x + y) mod W.
Definition sub x y := (x - y) mod W.
Definition bit (b : bool) : Z := if b then 1 else 0.

Lemma mod_up : forall s, 0 <= s < 2 * W -> s mod W = if s <? W then s else s - W.
Proof.
  intros s Hs; destruct (Z.ltb_spec s W).
  - apply Z.mod_small; lia.
  - symmetry; apply Z.mod_unique with (q := 1); unfold W in *; lia.
Qed.

Lemma mod_down : forall s, - W <= s < W -> s mod W = if s <? 0 then s + W else s.
Proof.
  intros s Hs; destruct (Z.ltb_spec s 0).
  - symmetry; apply Z.mod_unique with (q := -1); unfold W in *; lia.
  - apply Z.mod_small; lia.
Qed.

Ltac cases :=
  repeat match goal with
         | |- context [?x <? ?y] => destruct (Z.ltb_spec x y)
         | |- context [?x =? ?y] => destruct (Z.eqb_spec x y)
         end; cbn [bit orb andb negb] in *; unfold W in *; lia.

(* One step of an addition with carry, the carry taken by comparing with either addend: grid.cpp's
   Add compares with the term, whole.cpp's AddMagnitudes with the first. *)
Theorem add_step :
  forall a b c, word a -> word b -> 0 <= c <= 1 ->
  let partial := add a b in
  let total := add partial c in
  let carry := bit ((partial <? a) || (total <? partial)) in
  word total /\ carry = bit ((partial <? b) || (total <? partial)) /\ 0 <= carry <= 1 /\
  total + W * carry = a + b + c.
Proof.
  intros a b c Ha Hb Hc partial total carry; unfold word in *; subst partial total carry.
  unfold add; rewrite (mod_up (a + b)) by (unfold W in *; lia).
  destruct (Z.ltb_spec (a + b) W); rewrite mod_up by (unfold W in *; lia); cases.
Qed.

(* One step of a subtraction with borrow: whole.cpp's SubtractMagnitudes. *)
Theorem sub_step :
  forall a b c, word a -> word b -> 0 <= c <= 1 ->
  let partial := sub a b in
  let total := sub partial c in
  let borrow := bit ((a <? b) || (partial <? c)) in
  word total /\ 0 <= borrow <= 1 /\ total - W * borrow = a - b - c.
Proof.
  intros a b c Ha Hb Hc partial total borrow; unfold word in *; subst partial total borrow.
  unfold sub; rewrite (mod_down (a - b)) by (unfold W in *; lia).
  destruct (Z.ltb_spec (a - b) 0); rewrite mod_down by (unfold W in *; lia); cases.
Qed.

(* ~w + carry, the carry on while the words are zero: grid.cpp's Negated and Negate. *)
Definition complement x := W - 1 - x.

Theorem negate_step :
  forall a c, word a -> 0 <= c <= 1 ->
  let total := add (complement a) c in
  let carry := bit (negb (c =? 0) && (total =? 0)) in
  word total /\ 0 <= carry <= 1 /\ total + W * carry = complement a + c.
Proof.
  intros a c Ha Hc total carry; unfold word, complement in *; subst total carry.
  unfold add; rewrite mod_up by (unfold W in *; lia).
  destruct (Z.ltb_spec (W - 1 - a + c) W);
    destruct (Z.eqb_spec c 0); destruct (Z.eqb_spec (W - 1 - a + c) 0);
    destruct (Z.eqb_spec (W - 1 - a + c - W) 0); cbn [negb andb bit]; unfold W in *; lia.
Qed.

(* A product of two words: High and Low of grid.cpp's Multiply where the compiler has 128 bit
   integers, the product shifted right 64 and truncated. *)
Definition high_of p := p / W.
Definition low_of p := p mod W.

Lemma multiply_any :
  forall V a b, 1 < V -> 0 <= a < V -> 0 <= b < V ->
  0 <= (a * b) / V <= V - 2 /\ 0 <= (a * b) mod V < V /\ (a * b) / V * V + (a * b) mod V = a * b.
Proof.
  intros V a b HV Ha Hb.
  assert (Hab : a * b <= (V - 1) * (V - 1)) by (apply Z.mul_le_mono_nonneg; lia).
  assert (0 <= a * b) by (apply Z.mul_nonneg_nonneg; lia).
  assert ((V - 1) * (V - 1) < V * (V - 1)) by nia.
  pose proof (Z.div_mod (a * b) V ltac:(lia)).
  pose proof (Z.mod_pos_bound (a * b) V ltac:(lia)).
  assert (a * b / V < V - 1) by (apply Z.div_lt_upper_bound; lia).
  assert (0 <= a * b / V) by (apply Z.div_pos; lia).
  repeat split; try lia.
Qed.

Theorem multiply :
  forall a b, word a -> word b ->
  word (high_of (a * b)) /\ word (low_of (a * b)) /\
  high_of (a * b) * W + low_of (a * b) = a * b /\ high_of (a * b) <= W - 2.
Proof.
  intros a b Ha Hb; unfold word, high_of, low_of in *.
  destruct (multiply_any W a b) as [H1 [H2 H3]]; [unfold W; lia | lia | lia |].
  repeat split; lia.
Qed.

(* grid.hpp's MultiplyHalves, every compiler's: the 32 bit halves of each word, their four
   products, each below 2^64, and the middle column. x & 0xFFFFFFFF is x mod 2^32, x >> 32 is
   x / 2^32, middle << 32 is middle 2^32 mod 2^64. *)
Definition H : Z := 4294967296.

Definition halves_high_in V a b :=
  let lowLow := (a mod V) * (b mod V) in
  let lowHigh := (a mod V) * (b / V) in
  let highLow := (a / V) * (b mod V) in
  let highHigh := (a / V) * (b / V) in
  let middle := lowLow / V + lowHigh mod V + highLow mod V in
  highHigh + lowHigh / V + highLow / V + middle / V.

Definition halves_low_in V a b :=
  let lowLow := (a mod V) * (b mod V) in
  let lowHigh := (a mod V) * (b / V) in
  let highLow := (a / V) * (b mod V) in
  let middle := lowLow / V + lowHigh mod V + highLow mod V in
  (middle * V) mod (V * V) + lowLow mod V.

Lemma halves_any :
  forall V a b, 0 < V -> 0 <= a < V * V -> 0 <= b < V * V ->
  halves_high_in V a b = (a * b) / (V * V) /\ halves_low_in V a b = (a * b) mod (V * V).
Proof.
  intros V a b HV Ha Hb; unfold halves_high_in, halves_low_in; cbv zeta.
  pose proof (Z.div_mod a V ltac:(lia)) as Ea.
  pose proof (Z.div_mod b V ltac:(lia)) as Eb.
  pose proof (Z.mod_pos_bound a V HV) as Ba0.
  pose proof (Z.mod_pos_bound b V HV) as Bb0.
  assert (Ba1 : 0 <= a / V < V) by (split; [apply Z.div_pos | apply Z.div_lt_upper_bound]; lia).
  assert (Bb1 : 0 <= b / V < V) by (split; [apply Z.div_pos | apply Z.div_lt_upper_bound]; lia).
  set (a0 := a mod V) in *; set (a1 := a / V) in *;
  set (b0 := b mod V) in *; set (b1 := b / V) in *.
  clearbody a0 a1 b0 b1; subst a b.
  assert (Bll : 0 <= a0 * b0 < V * V) by nia.
  assert (Blh : 0 <= a0 * b1 < V * V) by nia.
  assert (Bhl : 0 <= a1 * b0 < V * V) by nia.
  set (ll := a0 * b0) in *; set (lh := a0 * b1) in *; set (hl := a1 * b0) in *.
  assert (Dll : ll = a0 * b0) by reflexivity.
  assert (Dlh : lh = a0 * b1) by reflexivity.
  assert (Dhl : hl = a1 * b0) by reflexivity.
  clearbody ll lh hl.
  pose proof (Z.div_mod ll V ltac:(lia)) as Ell.
  pose proof (Z.div_mod lh V ltac:(lia)) as Elh.
  pose proof (Z.div_mod hl V ltac:(lia)) as Ehl.
  pose proof (Z.mod_pos_bound ll V HV).
  pose proof (Z.mod_pos_bound lh V HV).
  pose proof (Z.mod_pos_bound hl V HV).
  assert (0 <= ll / V < V) by (split; [apply Z.div_pos | apply Z.div_lt_upper_bound]; lia).
  set (q3 := ll / V) in *; set (r3 := ll mod V) in *;
  set (q1 := lh / V) in *; set (r1 := lh mod V) in *;
  set (q2 := hl / V) in *; set (r2 := hl mod V) in *.
  clearbody q1 r1 q2 r2 q3 r3.
  assert (Bm : 0 <= q3 + r1 + r2) by lia.
  pose proof (Z.div_mod (q3 + r1 + r2) V ltac:(lia)) as Em.
  pose proof (Z.mod_pos_bound (q3 + r1 + r2) V HV).
  assert (Hm : ((q3 + r1 + r2) * V) mod (V * V) = ((q3 + r1 + r2) mod V) * V).
  { rewrite (Z.mul_comm (q3 + r1 + r2) V), Z.mul_mod_distr_l by lia; ring. }
  rewrite Hm.
  set (q4 := (q3 + r1 + r2) / V) in *; set (r4 := (q3 + r1 + r2) mod V) in *.
  clearbody q4 r4.
  assert (E : (V * a1 + a0) * (V * b1 + b0)
              = (a1 * b1 + q1 + q2 + q4) * (V * V) + (r4 * V + r3)).
  { assert (E0 : (V * a1 + a0) * (V * b1 + b0) = V * V * (a1 * b1) + V * (lh + hl) + ll)
      by (rewrite Dll, Dlh, Dhl; ring).
    rewrite E0, Elh, Ehl, Ell.
    replace r4 with (q3 + r1 + r2 - V * q4) by lia; ring. }
  assert (0 <= r4 * V + r3 < V * V) by nia.
  split.
  - apply Z.div_unique with (r := r4 * V + r3); [lia | rewrite E; ring].
  - apply Z.mod_unique with (q := a1 * b1 + q1 + q2 + q4); [lia | rewrite E; ring].
Qed.

Theorem multiply_halves :
  forall a b, word a -> word b ->
  halves_high_in H a b = high_of (a * b) /\ halves_low_in H a b = low_of (a * b).
Proof.
  intros a b Ha Hb; unfold word, high_of, low_of in *.
  replace W with (H * H) by reflexivity.
  apply halves_any; unfold H, W in *; lia.
Qed.


(* Whether a + b wraps, for words: the sum taken modulo 2^64 is below either addend exactly when
   it does. *)
Lemma add_wrap :
  forall a b, word a -> word b ->
  (add a b <? a) = (W <=? a + b) /\ (add a b <? b) = (W <=? a + b) /\
  add a b = if W <=? a + b then a + b - W else a + b.
Proof.
  intros a b Ha Hb; unfold word, add in *.
  rewrite mod_up by (unfold W in *; lia).
  destruct (Z.ltb_spec (a + b) W); destruct (Z.leb_spec W (a + b)); lia.
Qed.

Lemma W_gt : 1 < W.
Proof. unfold W; lia. Qed.

(* From here W is a number like any other to the tactics, its value used only through W_gt. *)
Opaque W.

(* ---------------------------------------------------------------------------------------------- *)
(* Magnitudes of many words, lowest first: whole.cpp's Whole, its words in use.                   *)

Fixpoint value (l : list Z) : Z :=
  match l with
  | [] => 0
  | w :: r => w + W * value r
  end.

Definition words l := Forall word l.

Fixpoint base (n : nat) : Z :=
  match n with
  | O => 1
  | S n => W * base n
  end.

Lemma base_pos : forall n, 1 <= base n.
Proof. induction n; cbn [base]; [lia|]. pose proof W_gt; nia. Qed.

Lemma base_le : forall n k, (n <= k)%nat -> base n <= base k.
Proof.
  intros n k Hnk; induction Hnk; [lia|]. cbn [base]. pose proof W_gt; pose proof (base_pos m); nia.
Qed.

Lemma base_add : forall n k, base (n + k) = base n * base k.
Proof. induction n; intros k; cbn [base Nat.add]; [ring|]. rewrite IHn; ring. Qed.

Lemma value_bounds : forall l, words l -> 0 <= value l < base (length l).
Proof.
  induction l as [|w r IH]; simpl; intros Hl; [lia|].
  inversion Hl as [|? ? Hw Hr]; subst; specialize (IH Hr); unfold word in Hw.
  assert (W * value r <= W * (base (length r) - 1)) by (apply Z.mul_le_mono_nonneg_l; lia).
  pose proof W_gt; nia.
Qed.

Lemma value_cons : forall l, value l = hd 0 l + W * value (tl l).
Proof. destruct l; simpl; lia. Qed.

Lemma words_hd : forall l, words l -> word (hd 0 l).
Proof. intros [|w r] H; simpl; [unfold word; pose proof W_gt; lia | now inversion H]. Qed.

Lemma words_tl : forall l, words l -> words (tl l).
Proof. intros [|w r] H; simpl; [constructor | now inversion H]. Qed.

Lemma value_app : forall l k, value (l ++ k) = value l + base (length l) * value k.
Proof.
  induction l as [|w r IH]; intros k; cbn [app value length base]; [ring|]. rewrite IH; ring.
Qed.

Lemma words_app : forall l k, words l -> words k -> words (l ++ k).
Proof. intros; now apply Forall_app. Qed.

(* AddMagnitudes' loop, over the longer's words, a missing word zero. *)
Fixpoint add_loop (n : nat) (xs ys : list Z) (carry : Z) : list Z * Z :=
  match n with
  | O => ([], carry)
  | S n =>
      let one := hd 0 xs in
      let two := hd 0 ys in
      let partial := add one two in
      let total := add partial carry in
      let rest := add_loop n (tl xs) (tl ys) (bit ((partial <? one) || (total <? partial))) in
      (total :: fst rest, snd rest)
  end.

Theorem add_loop_ok :
  forall n xs ys carry, words xs -> words ys -> (length xs <= n)%nat -> (length ys <= n)%nat ->
  0 <= carry <= 1 ->
  let r := add_loop n xs ys carry in
  words (fst r) /\ length (fst r) = n /\ 0 <= snd r <= 1 /\
  value (fst r) + base n * snd r = value xs + value ys + carry.
Proof.
  induction n as [|n IH]; intros xs ys carry Hx Hy Lx Ly Hc; simpl.
  - destruct xs, ys; cbn [length] in Lx, Ly; try lia.
    cbn [fst snd value base length]; repeat split;
        [apply Forall_nil | lia | lia | destruct carry; lia].
  - destruct (add_step (hd 0 xs) (hd 0 ys) carry (words_hd _ Hx) (words_hd _ Hy) Hc)
      as [Ht [_ [Hb E]]].
    destruct (IH (tl xs) (tl ys) _ (words_tl _ Hx) (words_tl _ Hy)
                 ltac:(destruct xs; simpl in *; lia) ltac:(destruct ys; simpl in *; lia) Hb)
      as [Hw [Hl [Hc' E']]].
    repeat split; try constructor; simpl; try lia; try assumption.
    rewrite (value_cons xs), (value_cons ys).
    replace (W * base n * snd _) with (W * (base n * snd (add_loop n (tl xs) (tl ys)
      (bit ((add (hd 0 xs) (hd 0 ys) <? hd 0 xs)
            || (add (add (hd 0 xs) (hd 0 ys)) carry <? add (hd 0 xs) (hd 0 ys))))))) by ring.
    nia.
Qed.

(* AddMagnitudes: the carry out, if any, a word more. *)
Definition add_magnitudes xs ys :=
  let r := add_loop (Nat.max (length xs) (length ys)) xs ys 0 in
  if snd r =? 0 then fst r else fst r ++ [snd r].

Theorem add_magnitudes_ok :
  forall xs ys, words xs -> words ys ->
  words (add_magnitudes xs ys) /\ value (add_magnitudes xs ys) = value xs + value ys.
Proof.
  intros xs ys Hx Hy; unfold add_magnitudes.
  destruct (add_loop_ok (Nat.max (length xs) (length ys)) xs ys 0 Hx Hy
              ltac:(lia) ltac:(lia) ltac:(lia)) as [Hw [Hl [Hc E]]].
  destruct (Z.eqb_spec (snd (add_loop (Nat.max (length xs) (length ys)) xs ys 0)) 0) as [Z0|Z1].
  - rewrite Z0 in E; split; [assumption | lia].
  - split.
    + apply words_app; [assumption | constructor; [unfold word in *; pose proof W_gt;
        lia | constructor]].
    + rewrite value_app, Hl; simpl; lia.
Qed.

(* SubtractMagnitudes' loop, over the larger's words. *)
Fixpoint sub_loop (larger smaller : list Z) (borrow : Z) : list Z * Z :=
  match larger with
  | [] => ([], borrow)
  | x :: rest =>
      let taken := hd 0 smaller in
      let partial := sub x taken in
      let total := sub partial borrow in
      let r := sub_loop rest (tl smaller) (bit ((x <? taken) || (partial <? borrow))) in
      (total :: fst r, snd r)
  end.

Theorem sub_loop_ok :
  forall larger smaller borrow, words larger -> words smaller ->
  (length smaller <= length larger)%nat -> 0 <= borrow <= 1 ->
  let r := sub_loop larger smaller borrow in
  words (fst r) /\ length (fst r) = length larger /\ 0 <= snd r <= 1 /\
  value (fst r) - base (length larger) * snd r = value larger - value smaller - borrow.
Proof.
  induction larger as [|x rest IH]; intros smaller borrow Hx Hy L Hb; simpl.
  - destruct smaller; cbn [length] in L; [|lia].
    cbn [fst snd value base length]; repeat split;
        [apply Forall_nil | lia | lia | destruct borrow; lia].
  - inversion Hx as [|? ? Hw Hr]; subst.
    destruct (sub_step x (hd 0 smaller) borrow Hw (words_hd _ Hy) Hb) as [Ht [Hb' E]].
    destruct (IH (tl smaller) _ Hr (words_tl _ Hy)
                 ltac:(destruct smaller; simpl in *; lia) Hb') as [Hw' [Hl [Hc' E']]].
    repeat split; try constructor; simpl; try lia; try assumption.
    rewrite (value_cons smaller).
    set (s := snd _) in *; nia.
Qed.

(* SubtractMagnitudes: where the first is at least the second, its loop borrows nothing out. *)
Theorem subtract_ok :
  forall larger smaller, words larger -> words smaller ->
  (length smaller <= length larger)%nat -> value smaller <= value larger ->
  words (fst (sub_loop larger smaller 0)) /\
  value (fst (sub_loop larger smaller 0)) = value larger - value smaller.
Proof.
  intros larger smaller Hx Hy L Hle.
  destruct (sub_loop_ok larger smaller 0 Hx Hy L ltac:(lia)) as [Hw [Hl [Hb E]]].
  pose proof (value_bounds _ Hw) as B; rewrite Hl in B.
  pose proof (value_bounds _ Hx).
  split; [assumption|].
  destruct (Z.eq_dec (snd (sub_loop larger smaller 0)) 0) as [Z0|Z1].
  - rewrite Z0 in E; lia.
  - exfalso; assert (snd (sub_loop larger smaller 0) = 1) by lia.
    pose proof (value_bounds _ Hy); nia.
Qed.

(* Words top first, as Compare and Trim walk them. *)
Fixpoint top (r : list Z) : Z :=
  match r with
  | [] => 0
  | w :: r' => w * base (length r') + top r'
  end.

Lemma value_rev : forall r, value (rev r) = top r.
Proof.
  induction r as [|w r IH]; [reflexivity|].
  cbn [rev top]; rewrite value_app, length_rev, IH; cbn [value]; ring.
Qed.

Lemma value_top : forall l, value l = top (rev l).
Proof. intros l; rewrite <- value_rev, rev_involutive; reflexivity. Qed.

Lemma words_rev : forall l, words l -> words (rev l).
Proof. intros l H; now apply Forall_rev. Qed.

Lemma top_bounds : forall r, words r -> 0 <= top r < base (length r).
Proof.
  intros r H; rewrite <- value_rev; rewrite <- (length_rev r); apply value_bounds, words_rev, H.
Qed.

(* Trim: the top words that are zero dropped. *)
Fixpoint drop_zeros (r : list Z) : list Z :=
  match r with
  | [] => []
  | w :: r' => if w =? 0 then drop_zeros r' else r
  end.

Definition trim l := rev (drop_zeros (rev l)).

(* Trimmed: no zero word on top, as Whole keeps its magnitude. *)
Definition trimmed l := hd 1 (rev l) <> 0.

Lemma drop_zeros_top : forall r, top (drop_zeros r) = top r.
Proof.
  induction r as [|w r IH]; [reflexivity|]; cbn [drop_zeros top].
  destruct (Z.eqb_spec w 0); [subst; rewrite IH; ring | reflexivity].
Qed.

Lemma drop_zeros_hd : forall r, hd 1 (drop_zeros r) <> 0.
Proof.
  induction r as [|w r IH]; cbn [drop_zeros hd]; [lia|].
  destruct (Z.eqb_spec w 0); [exact IH | cbn [hd]; exact n].
Qed.

Lemma drop_zeros_suffix : forall r, exists p, r = p ++ drop_zeros r.
Proof.
  induction r as [|w r [p IH]]; [exists []; reflexivity|]; cbn [drop_zeros].
  destruct (w =? 0); [exists (w :: p); cbn; f_equal; exact IH | exists []; reflexivity].
Qed.

Theorem trim_ok :
  forall l, words l -> words (trim l) /\ trimmed (trim l) /\ value (trim l) = value l /\
                       (length (trim l) <= length l)%nat.
Proof.
  intros l H; unfold trim, trimmed.
  destruct (drop_zeros_suffix (rev l)) as [p E].
  assert (Hw : words (rev l)) by now apply words_rev.
  rewrite E in Hw; apply Forall_app in Hw as [_ Hd].
  repeat split.
  - now apply words_rev.
  - rewrite rev_involutive; apply drop_zeros_hd.
  - rewrite value_rev, drop_zeros_top, <- value_top; reflexivity.
  - assert (length (rev l) = (length p + length (drop_zeros (rev l)))%nat)
      by (rewrite E at 1; apply length_app).
    rewrite length_rev in *; lia.
Qed.

Lemma trimmed_low :
  forall l, words l -> trimmed l -> l <> [] -> base (length l - 1) <= value l.
Proof.
  intros l H T N; unfold trimmed in T; rewrite value_top, <- (length_rev l).
  assert (Hr : words (rev l)) by now apply words_rev.
  destruct (rev l) as [|w r] eqn:R.
  - destruct l; [contradiction | cbn in R; destruct (rev l); discriminate].
  - cbn [hd] in T; cbn [top length]; inversion Hr as [|? ? Hw Hr']; subst.
    pose proof (top_bounds r Hr'); unfold word in Hw.
    replace (S (length r) - 1)%nat with (length r) by lia; pose proof (base_pos (length r)); nia.
Qed.

(* Compare: the longer trimmed magnitude is larger; of equal lengths, the first word down that
   differs decides. *)
Fixpoint compare_down (r s : list Z) : Z :=
  match r, s with
  | x :: r', y :: s' => if x =? y then compare_down r' s' else if y <? x then 1 else -1
  | _, _ => 0
  end.

Definition compare_magnitudes xs ys :=
  if negb (Nat.eqb (length xs) (length ys)) then (if Nat.ltb (length ys) (length xs) then 1 else -1)
  else compare_down (rev xs) (rev ys).

Lemma compare_down_ok :
  forall r s, words r -> words s -> length r = length s ->
  compare_down r s = Z.sgn (top r - top s).
Proof.
  induction r as [|x r IH]; intros [|y s] Hr Hs L; cbn [length] in L; try discriminate;
    [reflexivity|].
  injection L as L; inversion Hr as [|? ? Hx Hr']; inversion Hs as [|? ? Hy Hs']; subst.
  cbn [compare_down top length]; rewrite L.
  pose proof (top_bounds r Hr'); pose proof (top_bounds s Hs'); rewrite L in *.
  destruct (Z.eqb_spec x y).
  - subst; rewrite IH by (assumption || reflexivity); f_equal; ring.
  - set (B := base (length s)) in *.
    destruct (Z.ltb_spec y x).
    + symmetry; apply Z.sgn_pos; nia.
    + symmetry; apply Z.sgn_neg; nia.
Qed.

Theorem compare_ok :
  forall xs ys, words xs -> words ys -> trimmed xs -> trimmed ys ->
  compare_magnitudes xs ys = Z.sgn (value xs - value ys).
Proof.
  intros xs ys Hx Hy Tx Ty; unfold compare_magnitudes.
  pose proof (value_bounds xs Hx); pose proof (value_bounds ys Hy).
  destruct (Nat.eqb_spec (length xs) (length ys)) as [L|L]; cbn [negb].
  - rewrite compare_down_ok.
    + rewrite <- !value_top; reflexivity.
    + apply words_rev; assumption.
    + apply words_rev; assumption.
    + rewrite !length_rev; exact L.
  - destruct (Nat.ltb_spec (length ys) (length xs)).
    + assert (Nx : xs <> []) by (intros ->; cbn in *; lia).
      pose proof (trimmed_low xs Hx Tx Nx).
      pose proof (base_le (length ys) (length xs - 1) ltac:(lia)).
      symmetry; apply Z.sgn_pos; lia.
    + assert (Ny : ys <> []) by (intros ->; cbn in *; lia).
      pose proof (trimmed_low ys Hy Ty Ny).
      pose proof (base_le (length xs) (length ys - 1) ltac:(lia)).
      symmetry; apply Z.sgn_neg; lia.
Qed.

(* One step of Times' inner loop: the product of two words, plus the word already there, plus the
   carry, as a word and a carry that is again a word. *)
Lemma mac_step :
  forall x y w c, word x -> word y -> word w -> word c ->
  let low := add (low_of (x * y)) w in
  let total := add low c in
  let carry := add (add (high_of (x * y)) (bit (low <? w))) (bit (total <? low)) in
  word total /\ word carry /\ total + W * carry = x * y + w + c.
Proof.
  intros x y w c Hx Hy Hw Hc low total carry.
  destruct (multiply x y Hx Hy) as [Hh [Hl [E Hm]]].
  set (hi := high_of (x * y)) in *; set (lo := low_of (x * y)) in *.
  destruct (add_wrap lo w Hl Hw) as [_ [B1 S1]].
  fold low in B1, S1.
  assert (Hlow : word low) by (unfold low, add, word; apply Z.mod_pos_bound; pose proof W_gt; lia).
  destruct (add_wrap low c Hlow Hc) as [B2 [_ S2]].
  fold total in B2, S2.
  assert (Htotal : word total)
    by (unfold total, add, word; apply Z.mod_pos_bound; pose proof W_gt; lia).
  unfold word in *.
  set (b1 := bit (W <=? lo + w)); set (b2 := bit (W <=? low + c)).
  assert (Lb : low = lo + w - W * b1) by (unfold b1; destruct (Z.leb_spec W (lo + w)); cbn; lia).
  assert (Tb : total = low + c - W * b2) by (unfold b2; destruct (Z.leb_spec W (low + c)); cbn;
      lia).
  assert (0 <= b1 <= 1) by (unfold b1; destruct (W <=? lo + w); cbn; lia).
  assert (0 <= b2 <= 1) by (unfold b2; destruct (W <=? low + c); cbn; lia).
  assert (Hxy : x * y <= (W - 1) * (W - 1)) by (apply Z.mul_le_mono_nonneg; lia).
  assert (Sum : (hi + b1 + b2) * W + total = x * y + w + c) by lia.
  assert (hi + b1 + b2 < W) by nia.
  assert (A1 : add hi b1 = hi + b1).
  { destruct (add_wrap hi b1 ltac:(unfold word; lia) ltac:(unfold word; pose proof W_gt; lia))
      as [_ [_ S]]; rewrite S; destruct (Z.leb_spec W (hi + b1)); lia. }
  assert (A2 : add (hi + b1) b2 = hi + b1 + b2).
  { destruct (add_wrap (hi + b1) b2 ltac:(unfold word; lia) ltac:(unfold word; pose proof W_gt;
      lia))
      as [_ [_ S]]; rewrite S; destruct (Z.leb_spec W (hi + b1 + b2)); lia. }
  unfold carry; rewrite B1, B2; change (bit (W <=? lo + w)) with b1;
    change (bit (W <=? low + c)) with b2; rewrite A1, A2; repeat split; lia.
Qed.

(* Times' inner loop for one word of the first factor: each word of the second times it added to
   the product's word there, and the carry out written to the word past them. *)
Fixpoint mac (x : Z) (ys acc : list Z) (carry : Z) : list Z :=
  match ys with
  | [] => [carry]
  | y :: ys' =>
      let w := hd 0 acc in
      let low := add (low_of (x * y)) w in
      let total := add low carry in
      total :: mac x ys' (tl acc) (add (add (high_of (x * y)) (bit (low <? w))) (bit (total <?
          low)))
  end.

Lemma mac_ok :
  forall ys x acc c, words ys -> word x -> words acc -> (length ys <= length acc)%nat -> word c ->
  words (mac x ys acc c) /\ length (mac x ys acc c) = S (length ys) /\
  value (mac x ys acc c) = value (firstn (length ys) acc) + x * value ys + c.
Proof.
  induction ys as [|y ys IH]; intros x acc c Hy Hx Ha L Hc.
  - cbn [mac length firstn value]; repeat split.
    + constructor; [assumption | constructor].
    + ring.
  - destruct acc as [|w acc]; cbn [length] in L; [lia|].
    inversion Hy as [|? ? Hy1 Hys]; inversion Ha as [|? ? Hw Has]; subst.
    destruct (mac_step x y w c Hx Hy1 Hw Hc) as [Ht [Hc' E]].
    cbn [mac hd tl].
    destruct (IH x acc _ Hys Hx Has ltac:(lia) Hc') as [Hm [Lm Vm]].
    repeat split.
    + constructor; assumption.
    + cbn [length]; rewrite Lm; reflexivity.
    + cbn [value length firstn]; rewrite Vm; nia.
Qed.

Lemma words_firstn : forall k l, words l -> words (firstn k l).
Proof.
  intros k l H; rewrite <- (firstn_skipn k l) in H; apply Forall_app in H; apply H.
Qed.

Lemma words_skipn : forall k l, words l -> words (skipn k l).
Proof.
  intros k l H; rewrite <- (firstn_skipn k l) in H; apply Forall_app in H; apply H.
Qed.

Definition zeros l := Forall (fun w => w = 0) l.

Lemma zeros_value : forall l, zeros l -> value l = 0.
Proof.
  induction l as [|w r IH]; intros H; [reflexivity|].
  inversion H; subst; cbn [value]; rewrite IH by assumption; ring.
Qed.

Lemma zeros_skipn : forall k l, zeros l -> zeros (skipn k l).
Proof.
  intros k l H; rewrite <- (firstn_skipn k l) in H; apply Forall_app in H; apply H.
Qed.

Lemma value_split :
  forall k l, value l = value (firstn k l) + base k * value (skipn k l).
Proof.
  intros k l.
  destruct (Nat.le_gt_cases k (length l)).
  - rewrite <- (firstn_skipn k l) at 1; rewrite value_app, length_firstn.
    replace (Nat.min k (length l)) with k by lia; reflexivity.
  - rewrite firstn_all2, skipn_all2 by lia; cbn [value]; ring.
Qed.

(* Times' outer loop: row k adds the k-th word of the first factor times the second, W^k up. *)
Fixpoint rows (xs ys : list Z) (k : nat) (product : list Z) : list Z :=
  match xs with
  | [] => product
  | x :: xs' =>
      rows xs' ys (S k)
           (firstn k product ++ mac x ys (skipn k product) 0 ++ skipn (k + length ys + 1) product)
  end.

Lemma rows_ok :
  forall xs ys k product, words xs -> words ys -> words product ->
  length product = (k + length xs + length ys)%nat -> zeros (skipn (k + length ys) product) ->
  let r := rows xs ys k product in
  words r /\ length r = length product /\
  value r = value product + base k * value xs * value ys.
Proof.
  induction xs as [|x xs IH]; intros ys k product Hx Hy Hp L Z0.
  - cbn [rows value]; repeat split; [assumption | ring].
  - inversion Hx as [|? ? Hx1 Hxs]; subst; cbn [rows].
    set (m := length ys) in *.
    cbn [length] in L.
    assert (Hs : words (skipn k product)) by now apply words_skipn.
    assert (Ls : (m <= length (skipn k product))%nat) by (rewrite length_skipn; lia).
    destruct (mac_ok ys x (skipn k product) 0 Hy Hx1 Hs Ls ltac:(unfold word; pose proof W_gt; lia))
      as [Hm [Lm Vm]].
    set (M := mac x ys (skipn k product) 0) in *.
    set (P' := firstn k product ++ M ++ skipn (k + m + 1) product).
    assert (LP' : length P' = length product).
    { unfold P'; rewrite !length_app, length_firstn, Lm, length_skipn; lia. }
    assert (ZP' : zeros (skipn (S k + m) P')).
    { unfold P'; rewrite skipn_app, length_firstn.
      rewrite skipn_all2 by (rewrite length_firstn; lia).
      rewrite skipn_app, Lm.
      rewrite skipn_all2 by lia; cbn [app].
      assert (E0 : (S k + m - Nat.min k (length product) - S (length ys))%nat = 0%nat) by lia.
      rewrite E0.
      cbn [skipn]. replace (k + m + 1)%nat with (S (k + m)) by lia.
      replace (skipn (S (k + m)) product) with (skipn 1 (skipn (k + m) product))
        by (rewrite skipn_skipn; f_equal; lia).
      now apply zeros_skipn. }
    assert (HP' : words P').
    { unfold P'; apply words_app; [now apply words_firstn |].
      apply words_app; [assumption | now apply words_skipn]. }
    destruct (IH ys (S k) P' Hxs Hy HP' ltac:(rewrite LP'; lia) ZP') as [Hr [Lr Vr]].
    repeat split; [assumption | rewrite Lr; exact LP' |].
    rewrite Vr.
    assert (VP' : value P' = value product + base k * x * value ys).
    { unfold P'; rewrite !value_app, length_firstn, Lm.
      replace (Nat.min k (length product)) with k by lia.
      rewrite (zeros_value (skipn (k + m + 1) product))
        by (replace (k + m + 1)%nat with (1 + (k + m))%nat by lia;
            rewrite <- skipn_skipn; now apply zeros_skipn).
      rewrite Vm.
      rewrite (value_split k product).
      rewrite (value_split m (skipn k product)).
      rewrite skipn_skipn, (zeros_value (skipn (m + k) product))
        by (replace (m + k)%nat with (k + m)%nat by lia; assumption).
      unfold m; ring. }
    rewrite VP'; cbn [base value]; ring.
Qed.

Definition times xs ys := trim (rows xs ys 0 (repeat 0 (length xs + length ys))).

Lemma zeros_repeat : forall n, zeros (repeat 0 n).
Proof. induction n; constructor; auto. Qed.

Theorem times_ok :
  forall xs ys, words xs -> words ys ->
  words (times xs ys) /\ trimmed (times xs ys) /\ value (times xs ys) = value xs * value ys /\
  (length (times xs ys) <= length xs + length ys)%nat.
Proof.
  intros xs ys Hx Hy; unfold times.
  assert (Hz : words (repeat 0 (length xs + length ys))).
  { apply Forall_forall; intros w Hw; apply repeat_spec in Hw; subst; unfold word;
      pose proof W_gt; lia. }
  destruct (rows_ok xs ys 0 _ Hx Hy Hz ltac:(rewrite repeat_length; lia)
              (zeros_skipn _ _ (zeros_repeat _))) as [Hr [Lr Vr]].
  destruct (trim_ok _ Hr) as [Hw [T [V L]]].
  rewrite repeat_length in Lr.
  repeat split; try assumption; [| lia].
  rewrite V, Vr, zeros_value by apply zeros_repeat; cbn [base]; ring.
Qed.

(* ---------------------------------------------------------------------------------------------- *)
(* whole.cpp's signed whole numbers: a trimmed magnitude and a sign.                              *)

Record whole := Whole { mag : list Z; neg : bool }.

Definition signed (a : whole) := if neg a then - value (mag a) else value (mag a).
Definition good (a : whole) := words (mag a) /\ trimmed (mag a).

(* Plus: like signs add the magnitudes; unlike, the smaller is taken from the larger, which gives
   the sign. *)
Definition plus (a b : whole) (negate : bool) : whole :=
  let bneg := xorb (neg b) negate in
  if Bool.eqb (neg a) bneg then Whole (trim (add_magnitudes (mag a) (mag b))) (neg a)
  else if 0 <=? compare_magnitudes (mag a) (mag b)
       then Whole (trim (fst (sub_loop (mag a) (mag b) 0))) (neg a)
       else Whole (trim (fst (sub_loop (mag b) (mag a) 0))) bneg.

Lemma trimmed_length :
  forall xs ys, words xs -> words ys -> trimmed xs -> trimmed ys -> value ys <= value xs ->
  (length ys <= length xs)%nat.
Proof.
  intros xs ys Hx Hy Tx Ty Le.
  destruct (Nat.le_gt_cases (length ys) (length xs)) as [|G]; [assumption|].
  exfalso; assert (Ny : ys <> []) by (intros ->; cbn in G; lia).
  pose proof (trimmed_low ys Hy Ty Ny); pose proof (value_bounds xs Hx).
  pose proof (base_le (length xs) (length ys - 1) ltac:(lia)); lia.
Qed.

Theorem plus_ok :
  forall a b negate, good a -> good b ->
  good (plus a b negate) /\
  signed (plus a b negate) = signed a + (if negate then - signed b else signed b).
Proof.
  intros [xs na] [ys nb] negate [Hx Tx] [Hy Ty]; unfold plus, signed, good; cbn [mag neg] in *.
  pose proof (value_bounds xs Hx); pose proof (value_bounds ys Hy).
  destruct (Bool.eqb_spec na (xorb nb negate)) as [E|E].
  - destruct (add_magnitudes_ok xs ys Hx Hy) as [Hw V].
    destruct (trim_ok _ Hw) as [Hw' [T [V' _]]]; cbn [mag neg].
    repeat split; try assumption.
    rewrite V', V; destruct na, nb, negate; cbn in E; try discriminate; lia.
  - rewrite compare_ok by assumption.
    destruct (Z.leb_spec 0 (Z.sgn (value xs - value ys))) as [G|G].
    + assert (Le : value ys <= value xs) by (destruct (Z.sgn_spec (value xs - value ys)); lia).
      destruct (subtract_ok xs ys Hx Hy (trimmed_length xs ys Hx Hy Tx Ty Le) Le) as [Hw V].
      destruct (trim_ok _ Hw) as [Hw' [T [V' _]]]; cbn [mag neg].
      repeat split; try assumption.
      rewrite V', V; destruct na, nb, negate; cbn in E; try (exfalso; apply E; reflexivity); lia.
    + assert (Lt : value xs < value ys) by (destruct (Z.sgn_spec (value xs - value ys)); lia).
      destruct (subtract_ok ys xs Hy Hx (trimmed_length ys xs Hy Hx Ty Tx ltac:(lia)) ltac:(lia))
        as [Hw V].
      destruct (trim_ok _ Hw) as [Hw' [T [V' _]]]; cbn [mag neg].
      repeat split; try assumption.
      rewrite V', V; destruct na, nb, negate; cbn in E |- *; try (exfalso; apply E; reflexivity);
        lia.
Qed.

(* Times: the magnitudes' product, negative where the signs differ. *)
Definition times_whole (a b : whole) := Whole (times (mag a) (mag b)) (xorb (neg a) (neg b)).

Theorem times_whole_ok :
  forall a b, good a -> good b ->
  good (times_whole a b) /\ signed (times_whole a b) = signed a * signed b.
Proof.
  intros [xs na] [ys nb] [Hx _] [Hy _]; unfold times_whole, signed, good; cbn [mag neg] in *.
  destruct (times_ok xs ys Hx Hy) as [Hw [T [V _]]].
  repeat split; try assumption.
  rewrite V; destruct na, nb; cbn; ring.
Qed.

(* Sign: 0 for no words, else the sign flag's. *)
Definition sign_whole (a : whole) : Z :=
  match mag a with
  | [] => 0
  | _ => if neg a then -1 else 1
  end.

Theorem sign_whole_ok : forall a, good a -> sign_whole a = Z.sgn (signed a).
Proof.
  intros [xs na] [Hx Tx]; unfold sign_whole, signed; cbn [mag neg] in *.
  destruct xs as [|w r] eqn:X; [destruct na; reflexivity|].
  rewrite <- X in *.
  assert (Nx : xs <> []) by (rewrite X; discriminate).
  pose proof (trimmed_low xs Hx Tx Nx); pose proof (base_pos (length xs - 1)).
  destruct na; [symmetry; apply Z.sgn_neg | symmetry; apply Z.sgn_pos]; lia.
Qed.

(* OrientWhole after reading the coordinates in grid units: differences from the third corner,
   then their cross product. *)
Definition minus (a b : whole) := plus a b true.

Definition cross (fx fy sx sy : whole) := plus (times_whole fx sy) (times_whole fy sx) true.

Definition orient_whole (ax ay bx by' cx cy : whole) :=
  sign_whole (cross (minus ax cx) (minus ay cy) (minus bx cx) (minus by' cy)).

Lemma minus_ok :
  forall a b, good a -> good b -> good (minus a b) /\ signed (minus a b) = signed a - signed b.
Proof.
  intros a b Ha Hb; destruct (plus_ok a b true Ha Hb) as [G V].
  split; [exact G | unfold minus; rewrite V; ring].
Qed.

Lemma cross_ok :
  forall fx fy sx sy, good fx -> good fy -> good sx -> good sy ->
  good (cross fx fy sx sy) /\
  signed (cross fx fy sx sy) = signed fx * signed sy - signed fy * signed sx.
Proof.
  intros fx fy sx sy Hfx Hfy Hsx Hsy; unfold cross.
  destruct (times_whole_ok fx sy Hfx Hsy) as [G1 V1].
  destruct (times_whole_ok fy sx Hfy Hsx) as [G2 V2].
  destruct (plus_ok _ _ true G1 G2) as [G V]; split; [assumption | rewrite V, V1, V2; ring].
Qed.

Theorem orient_whole_ok :
  forall ax ay bx by' cx cy, good ax -> good ay -> good bx -> good by' -> good cx -> good cy ->
  orient_whole ax ay bx by' cx cy =
  Z.sgn ((signed ax - signed cx) * (signed by' - signed cy)
         - (signed ay - signed cy) * (signed bx - signed cx)).
Proof.
  intros ax ay bx by' cx cy Hax Hay Hbx Hby Hcx Hcy; unfold orient_whole.
  destruct (minus_ok ax cx Hax Hcx) as [G1 V1].
  destruct (minus_ok ay cy Hay Hcy) as [G2 V2].
  destruct (minus_ok bx cx Hbx Hcx) as [G3 V3].
  destruct (minus_ok by' cy Hby Hcy) as [G4 V4].
  destruct (cross_ok _ _ _ _ G1 G2 G3 G4) as [G V].
  rewrite sign_whole_ok, V, V1, V2, V3, V4 by assumption; reflexivity.
Qed.

(* CircleWhole after reading: each corner's difference from the point, its term the lift times the
   cross product of the other two, the three terms summed in turn. *)
Definition term (ox oy nx ny bx by' : whole) :=
  times_whole (plus (times_whole ox ox) (times_whole oy oy) false) (cross nx ny bx by').

Definition circle_whole (ax ay bx by' cx cy dx dy : whole) :=
  let ax' := minus ax dx in let ay' := minus ay dy in
  let bx' := minus bx dx in let by'' := minus by' dy in
  let cx' := minus cx dx in let cy' := minus cy dy in
  sign_whole (plus (plus (term ax' ay' bx' by'' cx' cy') (term bx' by'' cx' cy' ax' ay') false)
                   (term cx' cy' ax' ay' bx' by'') false).

Lemma term_ok :
  forall ox oy nx ny bx by', good ox -> good oy -> good nx -> good ny -> good bx -> good by' ->
  good (term ox oy nx ny bx by') /\
  signed (term ox oy nx ny bx by') =
  (signed ox * signed ox + signed oy * signed oy) * (signed nx * signed by' - signed ny * signed
      bx).
Proof.
  intros ox oy nx ny bx by' Hox Hoy Hnx Hny Hbx Hby; unfold term.
  destruct (times_whole_ok ox ox Hox Hox) as [G1 V1].
  destruct (times_whole_ok oy oy Hoy Hoy) as [G2 V2].
  destruct (plus_ok _ _ false G1 G2) as [G3 V3].
  destruct (cross_ok nx ny bx by' Hnx Hny Hbx Hby) as [G4 V4].
  destruct (times_whole_ok _ _ G3 G4) as [G V].
  split; [assumption | rewrite V, V3, V4, V1, V2; ring].
Qed.

Definition incircle (ax ay bx by' cx cy : Z) :=
  (ax * ax + ay * ay) * (bx * cy - by' * cx) + (bx * bx + by' * by') * (cx * ay - cy * ax)
  + (cx * cx + cy * cy) * (ax * by' - ay * bx).

Theorem circle_whole_ok :
  forall ax ay bx by' cx cy dx dy,
  good ax -> good ay -> good bx -> good by' -> good cx -> good cy -> good dx -> good dy ->
  circle_whole ax ay bx by' cx cy dx dy =
  Z.sgn (incircle (signed ax - signed dx) (signed ay - signed dy) (signed bx - signed dx)
                  (signed by' - signed dy) (signed cx - signed dx) (signed cy - signed dy)).
Proof.
  intros ax ay bx by' cx cy dx dy Hax Hay Hbx Hby Hcx Hcy Hdx Hdy; unfold circle_whole.
  destruct (minus_ok ax dx Hax Hdx) as [G1 V1].
  destruct (minus_ok ay dy Hay Hdy) as [G2 V2].
  destruct (minus_ok bx dx Hbx Hdx) as [G3 V3].
  destruct (minus_ok by' dy Hby Hdy) as [G4 V4].
  destruct (minus_ok cx dx Hcx Hdx) as [G5 V5].
  destruct (minus_ok cy dy Hcy Hdy) as [G6 V6].
  destruct (term_ok _ _ _ _ _ _ G1 G2 G3 G4 G5 G6) as [T1 W1].
  destruct (term_ok _ _ _ _ _ _ G3 G4 G5 G6 G1 G2) as [T2 W2].
  destruct (term_ok _ _ _ _ _ _ G5 G6 G1 G2 G3 G4) as [T3 W3].
  destruct (plus_ok _ _ false T1 T2) as [S1 U1].
  destruct (plus_ok _ _ false S1 T3) as [S2 U2].
  rewrite sign_whole_ok by assumption; f_equal.
  rewrite U2, U1, W1, W2, W3, V1, V2, V3, V4, V5, V6; unfold incircle; ring.
Qed.

(* ---------------------------------------------------------------------------------------------- *)
(* grid.cpp: two's complement in two words (Pair, and Product for magnitudes) and in four (Wide). *)

Definition HalfW : Z := 9223372036854775808.

Lemma W_half : W = 2 * HalfW.
Proof. rewrite W_pow; reflexivity. Qed.

Record pair := Pair { phigh : Z; plow : Z }.

Definition pvalue p := plow p + W * phigh p.
Definition pword p := word (phigh p) /\ word (plow p).

(* Multiply, either way it is computed. *)
Definition product a b := Pair (high_of (a * b)) (low_of (a * b)).

Lemma product_ok : forall a b,
    word a -> word b -> pword (product a b) /\ pvalue (product a b) = a * b.
Proof.
  intros a b Ha Hb; destruct (multiply a b Ha Hb) as [H1 [H2 [E _]]].
  unfold pword, pvalue, product; cbn [phigh plow]; repeat split; try apply H1; try apply H2; lia.
Qed.

(* (x >> 63) != 0 of a word: its top bit. *)
Definition top_bit x := negb (x / HalfW =? 0).

Lemma top_bit_ok : forall x, word x -> top_bit x = (HalfW <=? x).
Proof.
  intros x Hx; unfold top_bit, word in *; rewrite W_half in Hx.
  destruct (Z.leb_spec HalfW x).
  - assert (D : x / HalfW = 1).
    { symmetry; apply Z.div_unique with (r := x - HalfW); unfold HalfW in *; lia. }
    rewrite D; reflexivity.
  - rewrite Z.div_small by (unfold HalfW in *; lia); reflexivity.
Qed.

Definition below p := top_bit (phigh p).

(* Negate: ~High + (Low == 0), ~Low + 1. *)
Definition negate (p : pair) (n : bool) : pair :=
  if n then Pair (add (complement (phigh p)) (bit (plow p =? 0))) (add (complement (plow p)) 1)
  else p.

Definition W2 := W * W.

Lemma W2_pos : 0 < W2.
Proof. unfold W2; pose proof W_gt; nia. Qed.

Lemma pvalue_bounds : forall p, pword p -> 0 <= pvalue p < W2.
Proof. intros [h l] [Hh Hl]; unfold pvalue, W2, word in *; cbn [phigh plow] in *; nia. Qed.

Lemma negate_ok :
  forall p n, pword p -> pword (negate p n) /\
  pvalue (negate p n) = if n then (W2 - pvalue p) mod W2 else pvalue p.
Proof.
  intros [h l] n [Hh Hl]; destruct n; [|split; [split; assumption | reflexivity]].
  unfold negate, pvalue, pword, complement, add, word in *; cbn [phigh plow] in *.
  pose proof W_gt; unfold W2.
  destruct (Z.eqb_spec l 0) as [L0|L0]; cbn [bit].
  - subst l. rewrite (mod_up (W - 1 - h + 1)) by lia.
    rewrite (mod_up (W - 1 - 0 + 1)) by lia.
    destruct (Z.ltb_spec (W - 1 - h + 1) W); destruct (Z.ltb_spec (W - 1 - 0 + 1) W); try lia.
    + rewrite Z.mod_small by nia; repeat split; nia.
    + assert (h = 0) by lia; subst h.
      replace (W * W - (0 + W * 0)) with (W * W) by ring; rewrite Z.mod_same by nia.
      repeat split; lia.
  - rewrite (mod_up (W - 1 - h + 0)) by lia; rewrite (mod_up (W - 1 - l + 1)) by lia.
    destruct (Z.ltb_spec (W - 1 - h + 0) W); destruct (Z.ltb_spec (W - 1 - l + 1) W); try lia.
    rewrite Z.mod_small by nia; repeat split; nia.
Qed.

(* A two's complement value of two words, from its words. *)
Definition tc2 p := if below p then pvalue p - W2 else pvalue p.

Lemma below_ok : forall p, pword p -> below p = (HalfW * W <=? pvalue p).
Proof.
  intros [h l] [Hh Hl]; unfold below, pvalue; cbn [phigh plow] in *.
  rewrite top_bit_ok by assumption; unfold word in *.
  destruct (Z.leb_spec HalfW h); destruct (Z.leb_spec (HalfW * W) (l + W * h)); try reflexivity;
    pose proof W_gt; nia.
Qed.

Lemma tc2_of :
  forall p D, pword p -> pvalue p = D mod W2 -> - (HalfW * W) <= D < HalfW * W -> tc2 p = D.
Proof.
  intros p D Hp V B; unfold tc2; rewrite below_ok by assumption.
  assert (E : W2 = 2 * (HalfW * W)) by (unfold W2; rewrite W_half at 1; ring).
  destruct (Z.leb_spec 0 D).
  - rewrite Z.mod_small in V by lia; rewrite V.
    destruct (Z.leb_spec (HalfW * W) D); lia.
  - assert (V' : pvalue p = D + W2).
    { rewrite V; symmetry; apply Z.mod_unique with (q := -1); lia. }
    rewrite V'; destruct (Z.leb_spec (HalfW * W) (D + W2)); lia.
Qed.

Lemma tc2_bounds : forall p, pword p -> - (HalfW * W) <= tc2 p < HalfW * W.
Proof.
  intros p Hp; unfold tc2; rewrite below_ok by assumption.
  pose proof (pvalue_bounds p Hp).
  assert (E : W2 = 2 * (HalfW * W)) by (unfold W2; rewrite W_half at 1; ring).
  destruct (Z.leb_spec (HalfW * W) (pvalue p)); lia.
Qed.

Lemma tc2_mod : forall p, pword p -> pvalue p = tc2 p mod W2.
Proof.
  intros p Hp; pose proof (pvalue_bounds p Hp); unfold tc2.
  destruct (below p).
  - apply Z.mod_unique with (q := -1); lia.
  - symmetry; apply Z.mod_small; lia.
Qed.

(* Minus: the words' differences, the low word's borrow taken from the high. *)
Definition minus2 (e s : pair) : pair :=
  Pair (sub (sub (phigh e) (phigh s)) (bit (plow e <? plow s))) (sub (plow e) (plow s)).

Lemma mod_W2_of :
  forall v D, 0 <= v < W2 -> (v = D \/ v = D + W2 \/ v = D + 2 * W2) -> v = D mod W2.
Proof.
  intros v D Hv [E|[E|E]].
  - apply Z.mod_unique with (q := 0); lia.
  - apply Z.mod_unique with (q := -1); lia.
  - apply Z.mod_unique with (q := -2); lia.
Qed.

Lemma minus2_ok :
  forall e s, pword e -> pword s -> pword (minus2 e s) /\
  pvalue (minus2 e s) = (pvalue e - pvalue s) mod W2.
Proof.
  intros [eh el] [sh sl] [Heh Hel] [Hsh Hsl]; unfold minus2, pvalue, pword, sub, word in *;
    cbn [phigh plow] in *.
  pose proof W_gt; assert (E2 : W2 = W * W) by reflexivity.
  repeat (match goal with
          | |- context [?x mod W] => rewrite (mod_down x) by lia
          | |- context [?x <? ?y] => destruct (Z.ltb_spec x y)
          end; cbn beta iota delta [bit] in * ).
  all: repeat split; try lia.
  all: apply mod_W2_of; [nia | nia].
Qed.

Definition absolute p := negate p (below p).

Lemma absolute_ok : forall p, pword p -> pword (absolute p) /\ pvalue (absolute p) = Z.abs (tc2 p).
Proof.
  intros p Hp; unfold absolute; destruct (negate_ok p (below p) Hp) as [Hw V].
  split; [assumption|]; rewrite V; unfold tc2.
  pose proof (pvalue_bounds p Hp); pose proof (tc2_bounds p Hp); unfold tc2 in *.
  rewrite below_ok in * by assumption.
  assert (E : W2 = 2 * (HalfW * W)) by (unfold W2; rewrite W_half at 1; ring).
  destruct (Z.leb_spec (HalfW * W) (pvalue p)).
  - rewrite Z.mod_small by lia; lia.
  - lia.
Qed.

(* Wide: four words. grid.cpp's Add takes the carry by comparing with the term. *)
Fixpoint add_loop_term (n : nat) (xs ys : list Z) (carry : Z) : list Z * Z :=
  match n with
  | O => ([], carry)
  | S n =>
      let one := hd 0 xs in
      let two := hd 0 ys in
      let partial := add one two in
      let total := add partial carry in
      let rest := add_loop_term n (tl xs) (tl ys) (bit ((partial <? two) || (total <? partial))) in
      (total :: fst rest, snd rest)
  end.

Lemma add_loop_term_eq :
  forall n xs ys c,
      words xs -> words ys -> 0 <= c <= 1 -> add_loop_term n xs ys c = add_loop n xs ys c.
Proof.
  induction n as [|n IH]; intros xs ys c Hx Hy Hc; [reflexivity|]; cbn [add_loop_term add_loop].
  destruct (add_step (hd 0 xs) (hd 0 ys) c (words_hd _ Hx) (words_hd _ Hy) Hc) as [_ [E [B _]]].
  rewrite <- E, IH; [reflexivity | apply words_tl; assumption | apply words_tl;
      assumption | exact B].
Qed.

Definition B4 := base 4.

Definition add_wide s t := fst (add_loop_term 4 s t 0).

Definition wide4 w := words w /\ length w = 4%nat.

Lemma value_mod_of :
  forall v D k, 0 <= v < B4 -> v + B4 * k = D -> v = D mod B4.
Proof. intros v D k Hv E; apply Z.mod_unique with (q := k); lia. Qed.

Lemma add_wide_ok :
  forall s t, wide4 s -> wide4 t -> wide4 (add_wide s t) /\
  value (add_wide s t) = (value s + value t) mod B4.
Proof.
  intros s t [Hs Ls] [Ht Lt]; unfold add_wide; rewrite add_loop_term_eq by (try assumption; lia).
  destruct (add_loop_ok 4 s t 0 Hs Ht ltac:(lia) ltac:(lia) ltac:(lia)) as [Hw [L [C E]]].
  split; [split; assumption|].
  pose proof (value_bounds _ Hw) as B; rewrite L in B.
  apply value_mod_of with (k := snd (add_loop 4 s t 0)); [exact B | unfold B4; lia].
Qed.

(* Negated: ~w + carry, word by word, the carry 1 at first. *)
Fixpoint neg_loop (l : list Z) (carry : Z) : list Z * Z :=
  match l with
  | [] => ([], carry)
  | w :: r =>
      let total := add (complement w) carry in
      let rest := neg_loop r (bit (negb (carry =? 0) && (total =? 0))) in
      (total :: fst rest, snd rest)
  end.

Lemma neg_loop_ok :
  forall l c, words l -> 0 <= c <= 1 ->
  words (fst (neg_loop l c)) /\ length (fst (neg_loop l c)) = length l /\
  0 <= snd (neg_loop l c) <= 1 /\
  value (fst (neg_loop l c)) + base (length l) * snd (neg_loop l c)
  = base (length l) - 1 - value l + c.
Proof.
  induction l as [|w r IH]; intros c Hl Hc; cbn [neg_loop fst snd].
  - cbn [length value base]; repeat split; try apply Forall_nil; lia.
  - inversion Hl as [|? ? Hw Hr]; subst.
    destruct (negate_step w c Hw Hc) as [Ht [Hb E]].
    destruct (IH _ Hr Hb) as [Hw' [L [Hc' E']]].
    repeat split; [constructor; assumption | cbn [length]; rewrite L; reflexivity | lia | lia |].
    cbn [value length base]; unfold complement in *; nia.
Qed.

Definition negated l := fst (neg_loop l 1).

Lemma negated_ok :
  forall w, wide4 w -> wide4 (negated w) /\ value (negated w) = (B4 - value w) mod B4.
Proof.
  intros w [Hw L]; unfold negated.
  destruct (neg_loop_ok w 1 Hw ltac:(lia)) as [Hn [Ln [C E]]].
  pose proof (value_bounds _ Hn) as B; rewrite Ln, L in B; rewrite L in E.
  split; [split; [assumption | lia]|].
  apply value_mod_of with (k := snd (neg_loop w 1)); [exact B | unfold B4; lia].
Qed.

Definition signed_wide (m : pair) (n : bool) :=
  let w := [plow m; phigh m; 0; 0] in if n then negated w else w.

Definition negative_wide w := top_bit (nth 3 w 0).

Definition magnitude_wide w :=
  let w' := if negative_wide w then negated w else w in Pair (nth 1 w' 0) (nth 0 w' 0).

Definition tc4 w := if negative_wide w then value w - B4 else value w.

Lemma B4_eq : B4 = W * W * W * W.
Proof. unfold B4; cbn [base]; ring. Qed.

Lemma wide4_split :
  forall w, wide4 w -> exists a b c d, w = [a; b; c; d] /\ word a /\ word b /\ word c /\ word d.
Proof.
  intros w [Hw L].
  destruct w as [|a [|b [|c [|d [|e r]]]]]; cbn [length] in L; try discriminate.
  inversion Hw as [|? ? Ha I1]; inversion I1 as [|? ? Hb I2]; inversion I2 as [|? ? Hc I3];
    inversion I3 as [|? ? Hd _].
  exists a, b, c, d; exact (conj eq_refl (conj Ha (conj Hb (conj Hc Hd)))).
Qed.

Lemma negative_wide_ok : forall w, wide4 w -> negative_wide w = (HalfW * (W * W * W) <=? value w).
Proof.
  intros w Hw; destruct (wide4_split w Hw) as [a [b [c [d [-> [Ha [Hb [Hc Hd]]]]]]]].
  unfold negative_wide; cbn [nth value]; rewrite top_bit_ok by assumption.
  unfold word in *; pose proof W_gt.
  assert (a + W * (b + W * c) < W * W * W) by nia.
  destruct (Z.leb_spec HalfW d); destruct (Z.leb_spec (HalfW * (W * W * W))
    (a + W * (b + W * (c + W * (d + W * 0))))); try reflexivity; nia.
Qed.

Lemma tc4_of :
  forall w D, wide4 w -> value w = D mod B4 -> - (HalfW * (W * W * W)) <= D < HalfW * (W * W * W) ->
  tc4 w = D.
Proof.
  intros w D Hw V B; unfold tc4; rewrite negative_wide_ok by assumption.
  assert (E : B4 = 2 * (HalfW * (W * W * W))) by (rewrite B4_eq, W_half at 1; ring).
  destruct (Z.leb_spec 0 D).
  - rewrite Z.mod_small in V by lia; rewrite V.
    destruct (Z.leb_spec (HalfW * (W * W * W)) D); lia.
  - assert (V' : value w = D + B4).
    { rewrite V; symmetry; apply Z.mod_unique with (q := -1); lia. }
    rewrite V'; destruct (Z.leb_spec (HalfW * (W * W * W)) (D + B4)); lia.
Qed.

Lemma tc4_mod : forall w, wide4 w -> value w = tc4 w mod B4.
Proof.
  intros w Hw; pose proof (value_bounds w (proj1 Hw)) as B; rewrite (proj2 Hw) in B.
  fold B4 in B; unfold tc4.
  destruct (negative_wide w).
  - apply Z.mod_unique with (q := -1); lia.
  - symmetry; apply Z.mod_small; lia.
Qed.

(* SignOf of a wide integer: its top bit, else whether any word is set. *)
Definition sign_wide w :=
  if negative_wide w then -1 else if existsb (fun x => negb (x =? 0)) w then 1 else 0.

Lemma value_zero : forall l, words l -> value l = 0 -> zeros l.
Proof.
  induction l as [|w r IH]; intros H V; [constructor|].
  inversion H as [|? ? Hw Hr]; subst; cbn [value] in V.
  pose proof (value_bounds r Hr); unfold word in Hw; pose proof W_gt.
  assert (w = 0 /\ value r = 0) as [-> Vr] by nia.
  constructor; [reflexivity | apply IH; assumption].
Qed.

Lemma existsb_zeros :
  forall l, zeros l <-> existsb (fun x => negb (x =? 0)) l = false.
Proof.
  induction l as [|w r IH]; [split; [reflexivity | constructor]|].
  cbn [existsb]; rewrite orb_false_iff; split.
  - intros H; inversion H; subst; split; [reflexivity | apply IH; assumption].
  - intros [H1 H2]; constructor; [destruct (Z.eqb_spec w 0); [assumption | discriminate] |].
    apply IH; assumption.
Qed.

Lemma sign_wide_ok : forall w, wide4 w -> sign_wide w = Z.sgn (tc4 w).
Proof.
  intros w Hw; unfold sign_wide, tc4.
  pose proof (value_bounds w (proj1 Hw)) as B; rewrite (proj2 Hw) in B.
  rewrite negative_wide_ok in * by assumption.
  assert (E : base 4 = 2 * (HalfW * (W * W * W))) by (cbn [base]; rewrite W_half at 1; ring).
  destruct (Z.leb_spec (HalfW * (W * W * W)) (value w)).
  - symmetry; apply Z.sgn_neg; fold B4 in B; unfold B4 in *; lia.
  - destruct (existsb (fun x => negb (x =? 0)) w) eqn:X.
    + symmetry; apply Z.sgn_pos.
      assert (value w <> 0).
      { intros V0; apply value_zero in V0; [| exact (proj1 Hw)].
        apply existsb_zeros in V0; congruence. }
      lia.
    + rewrite (zeros_value w (proj2 (existsb_zeros w) X)); reflexivity.
Qed.

Lemma low_two :
  forall w, wide4 w -> value w < W * W -> pword (Pair (nth 1 w 0) (nth 0 w 0)) /\
  pvalue (Pair (nth 1 w 0) (nth 0 w 0)) = value w.
Proof.
  intros w Hw V; destruct (wide4_split w Hw) as [a [b [c [d [-> [Ha [Hb [Hc Hd]]]]]]]].
  unfold pword, pvalue, word in *; cbn [nth phigh plow value] in *; pose proof W_gt.
  replace (a + W * (b + W * (c + W * (d + W * 0)))) with ((a + W * b) + (W * W) * (c + W * d))
    in * by ring.
  assert (0 <= W * b) by nia; assert (0 <= W * d) by nia.
  set (X := W * W) in *; assert (1 < X) by (unfold X; nia).
  assert (c + W * d < 1) by nia.
  assert (c = 0 /\ d = 0) as [-> ->] by nia.
  repeat split; try lia.
Qed.

Lemma magnitude_wide_ok :
  forall w, wide4 w -> Z.abs (tc4 w) < W * W ->
  pword (magnitude_wide w) /\ pvalue (magnitude_wide w) = Z.abs (tc4 w).
Proof.
  intros w Hw A; unfold magnitude_wide, tc4 in *.
  pose proof (value_bounds w (proj1 Hw)) as B; rewrite (proj2 Hw) in B; fold B4 in B.
  destruct (negative_wide w) eqn:N.
  - destruct (negated_ok w Hw) as [Hn V].
    pose proof N as N'; rewrite negative_wide_ok in N' by assumption; apply Z.leb_le in N'.
    assert (0 < HalfW * (W * W * W)) by (pose proof W_gt; unfold HalfW; nia).
    assert (V' : value (negated w) = B4 - value w) by (rewrite V; apply Z.mod_small; lia).
    destruct (low_two (negated w) Hn ltac:(lia)) as [P E]; split; [exact P | lia].
  - destruct (low_two w Hw ltac:(lia)) as [P E]; split; [exact P | lia].
Qed.

(* grid.cpp's Times of two magnitudes of two words: four products, two additions. *)
Definition times_wide (a b : pair) :=
  let ll := product (plow a) (plow b) in
  let lh := product (plow a) (phigh b) in
  let hl := product (phigh a) (plow b) in
  let hh := product (phigh a) (phigh b) in
  add_wide (add_wide [plow ll; phigh ll; plow hh; phigh hh] [0; plow lh; phigh lh; 0])
           [0; plow hl; phigh hl; 0].

Lemma wide4_of :
  forall a b c d, word a -> word b -> word c -> word d -> wide4 [a; b; c; d].
Proof.
  intros; split; [repeat (apply Forall_cons; [assumption|]); apply Forall_nil | reflexivity].
Qed.

Lemma word_zero : word 0.
Proof. unfold word; pose proof W_gt; lia. Qed.

Theorem times_wide_ok :
  forall a b,
      pword a -> pword b -> wide4 (times_wide a b) /\ value (times_wide a b) = pvalue a * pvalue b.
Proof.
  intros [ah al] [bh bl] [Hah Hal] [Hbh Hbl]; unfold times_wide, pvalue; cbn [phigh plow] in *.
  destruct (product_ok al bl Hal Hbl) as [[P1 P2] V1].
  destruct (product_ok al bh Hal Hbh) as [[P3 P4] V2].
  destruct (product_ok ah bl Hah Hbl) as [[P5 P6] V3].
  destruct (product_ok ah bh Hah Hbh) as [[P7 P8] V4].
  unfold pvalue in *.
  set (ll := product al bl) in *; set (lh := product al bh) in *;
  set (hl := product ah bl) in *; set (hh := product ah bh) in *.
  pose proof word_zero as Z0.
  destruct (add_wide_ok [plow ll; phigh ll; plow hh; phigh hh] [0; plow lh; phigh lh; 0]
              (wide4_of _ _ _ _ P2 P1 P8 P7) (wide4_of _ _ _ _ Z0 P4 P3 Z0)) as [S1 U1].
  destruct (add_wide_ok _ [0; plow hl; phigh hl; 0] S1 (wide4_of _ _ _ _ Z0 P6 P5 Z0))
    as [S2 U2].
  split; [exact S2|].
  rewrite U2, U1, Zplus_mod_idemp_l; cbn [value].
  unfold word in *; pose proof W_gt.
  assert (0 <= (al + W * ah) * (bl + W * bh) < B4).
  { rewrite B4_eq; split; [nia|].
    assert (al + W * ah < W * W) by nia; assert (bl + W * bh < W * W) by nia; nia. }
  rewrite Z.mod_small; [nia|].
  split; [nia|].
  replace (plow ll + W * (phigh ll + W * (plow hh + W * (phigh hh + W * 0))) +
           (0 + W * (plow lh + W * (phigh lh + W * (0 + W * 0)))) +
           (0 + W * (plow hl + W * (phigh hl + W * (0 + W * 0)))))
    with ((al + W * ah) * (bl + W * bh)) by nia.
  lia.
Qed.

Lemma below_tc2 : forall p, pword p -> below p = (tc2 p <? 0).
Proof.
  intros p Hp; unfold tc2; pose proof (pvalue_bounds p Hp).
  rewrite below_ok by assumption.
  assert (E : W2 = 2 * (HalfW * W)) by (unfold W2; rewrite W_half at 1; ring).
  destruct (Z.leb_spec (HalfW * W) (pvalue p)); [symmetry; apply Z.ltb_lt; lia|].
  symmetry; apply Z.ltb_ge; lia.
Qed.

(* Signed: a magnitude of two words, taken below zero or not, in four. *)
Lemma signed_wide_ok :
  forall m n, pword m -> wide4 (signed_wide m n) /\
  value (signed_wide m n) = (if n then - pvalue m else pvalue m) mod B4.
Proof.
  intros [h l] n [Hh Hl]; unfold signed_wide, pvalue; cbn [phigh plow] in *.
  pose proof word_zero as Z0.
  assert (Hw : wide4 [l; h; 0; 0]) by (apply wide4_of; assumption).
  assert (V : value [l; h; 0; 0] = l + W * h) by (cbn [value]; ring).
  assert (B : 0 <= l + W * h < B4).
  { rewrite B4_eq; unfold word in *; pose proof W_gt.
    assert (l + W * h < W * W) by nia; assert (W * W <= W * W * W * W) by nia; nia. }
  destruct n.
  - destruct (negated_ok _ Hw) as [Hn Vn]; split; [assumption|].
    rewrite Vn, V.
    destruct (Z.eq_dec (l + W * h) 0) as [E|E].
    + rewrite E, Z.sub_0_r, Z.mod_same by lia; reflexivity.
    + rewrite Z.mod_small by lia; apply Z.mod_unique with (q := -1); lia.
  - split; [assumption|]; rewrite V; symmetry; apply Z.mod_small; lia.
Qed.

(* Term's word(value): an int64 read as its magnitude in a word, value < 0 ? 0 - (uint64)value
   : (uint64)value, the conversion to std::uint64_t modulo 2^64. *)
Definition word_of (v : Z) := if v <? 0 then sub 0 (v mod W) else v mod W.

Lemma word_of_abs : forall v, Z.abs v < W -> word_of v = Z.abs v /\ word (word_of v).
Proof.
  intros v A; unfold word_of, sub, word; pose proof W_gt.
  destruct (Z.ltb_spec v 0).
  - assert (E1 : v mod W = v + W) by (symmetry; apply Z.mod_unique with (q := -1); lia).
    rewrite E1.
    assert (E2 : (0 - (v + W)) mod W = - v) by (symmetry; apply Z.mod_unique with (q := -1); lia).
    rewrite E2; lia.
  - rewrite Z.mod_small by lia; lia.
Qed.

(* A product of two int64 below 2^64 in magnitude, as Term forms it: the magnitudes' product,
   below zero where their signs differ. *)
Lemma signed_product :
  forall a b n, Z.abs a < W -> Z.abs b < W -> n = xorb (a <? 0) (b <? 0) ->
  wide4 (signed_wide (product (word_of a) (word_of b)) n) /\
  value (signed_wide (product (word_of a) (word_of b)) n) = (a * b) mod B4.
Proof.
  intros a b n Aa Ab N.
  destruct (word_of_abs a Aa) as [La Wa]; destruct (word_of_abs b Ab) as [Lb Wb].
  destruct (product_ok _ _ Wa Wb) as [Hp Vp].
  destruct (signed_wide_ok _ n Hp) as [Hs Vs]; split; [assumption|].
  rewrite Vs, Vp, La, Lb, N.
  destruct (Z.ltb_spec a 0); destruct (Z.ltb_spec b 0); cbn [xorb]; f_equal.
  - rewrite Z.abs_neq, (Z.abs_neq b) by lia; ring.
  - rewrite Z.abs_neq, (Z.abs_eq b) by lia; ring.
  - rewrite Z.abs_eq, (Z.abs_neq b) by lia; ring.
  - rewrite Z.abs_eq, (Z.abs_eq b) by lia; ring.
Qed.

(* Adding in four words is adding modulo 2^256. *)
Lemma add_wide_mod :
  forall s t S T, wide4 s -> wide4 t -> value s = S mod B4 -> value t = T mod B4 ->
  wide4 (add_wide s t) /\ value (add_wide s t) = (S + T) mod B4.
Proof.
  intros s t S T Hs Ht Vs Vt; destruct (add_wide_ok s t Hs Ht) as [H V]; split; [assumption|].
  rewrite V, Vs, Vt, <- Zplus_mod; reflexivity.
Qed.

Lemma half_bound : HalfW * (W * W * W) = 2 ^ 255.
Proof. rewrite W_pow; reflexivity. Qed.

Lemma W2_bound : W * W = 2 ^ 128.
Proof. rewrite W_pow; reflexivity. Qed.

(* grid.cpp's Term: a corner's lift times the minor of the two after it, each difference an
   int64 below 2^62 in magnitude. *)
Definition term_grid (ox oy nx ny bx by' : Z) :=
  let lift := add_wide (signed_wide (product (word_of ox) (word_of ox)) false)
                       (signed_wide (product (word_of oy) (word_of oy)) false) in
  let minor := add_wide (signed_wide (product (word_of nx) (word_of by'))
                                     (xorb (nx <? 0) (by' <? 0)))
                        (signed_wide (product (word_of ny) (word_of bx))
                                     (Bool.eqb (ny <? 0) (bx <? 0))) in
  let whole := times_wide (magnitude_wide lift) (magnitude_wide minor) in
  if negative_wide minor then negated whole else whole.

Lemma small_W : forall v, Z.abs v < 2 ^ 62 -> Z.abs v < W.
Proof. intros v A; rewrite W_pow; lia. Qed.

Theorem term_grid_ok :
  forall ox oy nx ny bx by',
  Z.abs ox < 2 ^ 62 -> Z.abs oy < 2 ^ 62 -> Z.abs nx < 2 ^ 62 -> Z.abs ny < 2 ^ 62 ->
  Z.abs bx < 2 ^ 62 -> Z.abs by' < 2 ^ 62 ->
  wide4 (term_grid ox oy nx ny bx by') /\
  value (term_grid ox oy nx ny bx by') = ((ox * ox + oy * oy) * (nx * by' - ny * bx)) mod B4.
Proof.
  intros ox oy nx ny bx by' Aox Aoy Anx Any Abx Aby; unfold term_grid; cbv zeta.
  pose proof (small_W _ Aox) as Wox; pose proof (small_W _ Aoy) as Woy;
  pose proof (small_W _ Anx) as Wnx; pose proof (small_W _ Any) as Wny;
  pose proof (small_W _ Abx) as Wbx; pose proof (small_W _ Aby) as Wby.
  assert (Xo : false = xorb (ox <? 0) (ox <? 0)) by (destruct (ox <? 0); reflexivity).
  assert (Yo : false = xorb (oy <? 0) (oy <? 0)) by (destruct (oy <? 0); reflexivity).
  destruct (signed_product ox ox false Wox Wox Xo) as [L1 V1].
  destruct (signed_product oy oy false Woy Woy Yo) as [L2 V2].
  destruct (add_wide_mod _ _ _ _ L1 L2 V1 V2) as [Lw Lv].
  destruct (signed_product nx by' _ Wnx Wby eq_refl) as [M1 U1].
  assert (M2 : wide4 (signed_wide (product (word_of ny) (word_of bx))
                                  (Bool.eqb (ny <? 0) (bx <? 0)))
               /\ value (signed_wide (product (word_of ny) (word_of bx))
                                     (Bool.eqb (ny <? 0) (bx <? 0))) = (- (ny * bx)) mod B4).
  { destruct (word_of_abs ny Wny) as [La Wa]; destruct (word_of_abs bx Wbx) as [Lb Wb].
    destruct (product_ok _ _ Wa Wb) as [Hp Vp].
    destruct (signed_wide_ok _ (Bool.eqb (ny <? 0) (bx <? 0)) Hp) as [Hs Vs].
    split; [assumption|]; rewrite Vs, Vp, La, Lb.
    destruct (Z.ltb_spec ny 0); destruct (Z.ltb_spec bx 0); cbn [Bool.eqb]; f_equal.
    - rewrite Z.abs_neq, (Z.abs_neq bx) by lia; ring.
    - rewrite Z.abs_neq, (Z.abs_eq bx) by lia; ring.
    - rewrite Z.abs_eq, (Z.abs_neq bx) by lia; ring.
    - rewrite Z.abs_eq, (Z.abs_eq bx) by lia; ring. }
  destruct M2 as [M2 U2].
  destruct (add_wide_mod _ _ _ _ M1 M2 U1 U2) as [Mw Mv].
  set (lift := add_wide _ _) in *; set (minor := add_wide _ (signed_wide _ (Bool.eqb _ _))) in *.
  set (L := ox * ox + oy * oy) in *.
  set (M := nx * by' + - (ny * bx)) in *.
  assert (BL : 0 <= L < 2 ^ 125) by (unfold L; nia).
  assert (BM : Z.abs M < 2 ^ 125) by (unfold M; nia).
  assert (TL : tc4 lift = L) by (apply tc4_of; [assumption | assumption | rewrite half_bound; lia]).
  assert (TM : tc4 minor = M)
    by (apply tc4_of; [assumption | assumption | rewrite half_bound; lia]).
  destruct (magnitude_wide_ok lift Lw ltac:(rewrite TL, W2_bound; lia)) as [PL VL].
  destruct (magnitude_wide_ok minor Mw ltac:(rewrite TM, W2_bound; lia)) as [PM VM].
  destruct (times_wide_ok _ _ PL PM) as [Tw Tv].
  rewrite VL, VM, TL, TM in Tv.
  assert (NM : negative_wide minor = (M <? 0)).
  { unfold tc4 in TM; rewrite negative_wide_ok in * by assumption.
    pose proof (value_bounds minor (proj1 Mw)) as B; rewrite (proj2 Mw) in B; fold B4 in B.
    assert (E : B4 = 2 * (HalfW * (W * W * W))) by (rewrite B4_eq, W_half at 1; ring).
    destruct (Z.leb_spec (HalfW * (W * W * W)) (value minor)); symmetry;
      [apply Z.ltb_lt | apply Z.ltb_ge]; lia. }
  rewrite NM; replace (nx * by' - ny * bx) with M by (unfold M; ring).
  destruct (Z.ltb_spec M 0).
  - destruct (negated_ok _ Tw) as [Nw Nv]; split; [assumption|].
    rewrite Nv, Tv, Z.abs_eq, Z.abs_neq by lia.
    replace (B4 - L * - M) with (L * M + 1 * B4) by ring.
    rewrite Z.mod_add by (rewrite B4_eq; pose proof W_gt; nia); reflexivity.
  - split; [assumption|]; rewrite Tv, !Z.abs_eq by lia.
    symmetry; apply Z.mod_small; rewrite B4_eq, W_pow; nia.
Qed.

(* Differences of grid units: Minus, exact while they stay below 2^127 in magnitude. *)
Lemma minus2_tc :
  forall e s, pword e -> pword s -> Z.abs (tc2 e - tc2 s) < HalfW * W ->
  pword (minus2 e s) /\ tc2 (minus2 e s) = tc2 e - tc2 s.
Proof.
  intros e s He Hs A; destruct (minus2_ok e s He Hs) as [Hm V]; split; [assumption|].
  apply tc2_of; [assumption | | lia].
  rewrite V, (tc2_mod e He), (tc2_mod s Hs), <- Zminus_mod; reflexivity.
Qed.

Lemma square_bound : forall a X, Z.abs a < X -> 0 <= a * a < X * X.
Proof. intros a X A; destruct (Z.abs_spec a) as [[_ E]|[_ E]]; rewrite E in A; nia. Qed.

Lemma product_bound : forall a b X, Z.abs a < X -> Z.abs b < X -> Z.abs (a * b) < X * X.
Proof.
  intros a b X A B; rewrite Z.abs_mul.
  pose proof (Z.abs_nonneg a); pose proof (Z.abs_nonneg b); nia.
Qed.

Lemma term_bound :
  forall a b c d e f, Z.abs a < 2 ^ 62 -> Z.abs b < 2 ^ 62 -> Z.abs c < 2 ^ 62 ->
  Z.abs d < 2 ^ 62 -> Z.abs e < 2 ^ 62 -> Z.abs f < 2 ^ 62 ->
  Z.abs ((a * a + b * b) * (c * f - d * e)) < 2 ^ 250.
Proof.
  intros a b c d e f Ha Hb Hc Hd He Hf.
  pose proof (square_bound a _ Ha); pose proof (square_bound b _ Hb).
  pose proof (product_bound c f _ Hc Hf); pose proof (product_bound d e _ Hd He).
  set (X := 2 ^ 62) in *.
  assert (Z.abs (c * f - d * e) < 2 * (X * X)) by lia.
  rewrite Z.abs_mul, (Z.abs_eq (a * a + b * b)) by lia.
  replace (2 ^ 250) with ((2 * (X * X)) * (2 * (X * X))) by (unfold X; reflexivity).
  pose proof (Z.abs_nonneg (c * f - d * e)); nia.
Qed.

(* OnGrid after reading: each corner's difference from the point, in int64, exact below 2^62,
   the three terms summed. *)
Definition on_grid (ax ay bx by' cx cy dx dy : Z) :=
  let ax' := ax - dx in let ay' := ay - dy in
  let bx' := bx - dx in let by'' := by' - dy in
  let cx' := cx - dx in let cy' := cy - dy in
  sign_wide (add_wide (add_wide (term_grid ax' ay' bx' by'' cx' cy')
                                (term_grid bx' by'' cx' cy' ax' ay'))
                      (term_grid cx' cy' ax' ay' bx' by'')).

Theorem on_grid_ok :
  forall ax ay bx by' cx cy dx dy,
  Z.abs ax < 2 ^ 61 -> Z.abs ay < 2 ^ 61 -> Z.abs bx < 2 ^ 61 -> Z.abs by' < 2 ^ 61 ->
  Z.abs cx < 2 ^ 61 -> Z.abs cy < 2 ^ 61 -> Z.abs dx < 2 ^ 61 -> Z.abs dy < 2 ^ 61 ->
  on_grid ax ay bx by' cx cy dx dy =
  Z.sgn (incircle (ax - dx) (ay - dy) (bx - dx) (by' - dy) (cx - dx) (cy - dy)).
Proof.
  intros ax ay bx by' cx cy dx dy Hax Hay Hbx Hby Hcx Hcy Hdx Hdy; unfold on_grid; cbv zeta.
  assert (A1 : Z.abs (ax - dx) < 2 ^ 62) by lia; assert (A2 : Z.abs (ay - dy) < 2 ^ 62) by lia.
  assert (A3 : Z.abs (bx - dx) < 2 ^ 62) by lia; assert (A4 : Z.abs (by' - dy) < 2 ^ 62) by lia.
  assert (A5 : Z.abs (cx - dx) < 2 ^ 62) by lia; assert (A6 : Z.abs (cy - dy) < 2 ^ 62) by lia.
  generalize dependent (cy - dy); generalize dependent (cx - dx);
  generalize dependent (by' - dy); generalize dependent (bx - dx);
  generalize dependent (ay - dy); generalize dependent (ax - dx).
  intros a1 A1 a2 A2 b1 A3 b2 A4 c1 A5 c2 A6.
  destruct (term_grid_ok _ _ _ _ _ _ A1 A2 A3 A4 A5 A6) as [W1 V1].
  destruct (term_grid_ok _ _ _ _ _ _ A3 A4 A5 A6 A1 A2) as [W2 V2].
  destruct (term_grid_ok _ _ _ _ _ _ A5 A6 A1 A2 A3 A4) as [W3 V3].
  destruct (add_wide_mod _ _ _ _ W1 W2 V1 V2) as [W12 V12].
  destruct (add_wide_mod _ _ _ _ W12 W3 V12 V3) as [W123 V123].
  rewrite sign_wide_ok by assumption; f_equal.
  apply tc4_of; [assumption | rewrite V123; f_equal; unfold incircle; ring |].
  rewrite half_bound.
  pose proof (term_bound _ _ _ _ _ _ A1 A2 A3 A4 A5 A6) as B1.
  pose proof (term_bound _ _ _ _ _ _ A3 A4 A5 A6 A1 A2) as B2.
  pose proof (term_bound _ _ _ _ _ _ A5 A6 A1 A2 A3 A4) as B3.
  clear - B1 B2 B3; unfold incircle.
  generalize dependent ((c1 * c1 + c2 * c2) * (a1 * b2 - a2 * b1)).
  generalize dependent ((b1 * b1 + b2 * b2) * (c1 * a2 - c2 * a1)).
  generalize dependent ((a1 * a1 + a2 * a2) * (b1 * c2 - b2 * c1)).
  intros t1 B1 t2 B2 t3 B3; lia.
Qed.

(* OrientGrid's sign of one product from its magnitude and its factors' signs. *)
Definition sign_product (m : list Z) (f s : pair) :=
  if forallb (fun x => x =? 0) m then 0 else if xorb (below f) (below s) then -1 else 1.

Lemma forallb_zeros : forall l, forallb (fun x => x =? 0) l = true <-> zeros l.
Proof.
  induction l as [|w r IH]; [split; [constructor | reflexivity]|].
  cbn [forallb]; rewrite andb_true_iff, IH; split.
  - intros [A B]; constructor; [apply Z.eqb_eq; assumption | assumption].
  - intros H; inversion H; subst; split; [apply Z.eqb_refl | assumption].
Qed.

Definition orient_grid (ax ay bx by' cx cy : pair) :=
  let fx := minus2 ax cx in let fy := minus2 ay cy in
  let sx := minus2 bx cx in let sy := minus2 by' cy in
  let left := times_wide (absolute fx) (absolute sy) in
  let right := times_wide (absolute fy) (absolute sx) in
  let ls := sign_product left fx sy in
  let rs := sign_product right fy sx in
  if negb (ls =? rs) then (if rs <? ls then 1 else -1) else ls * compare_down (rev left) (rev
      right).

Definition unit125 p := pword p /\ Z.abs (tc2 p) < 2 ^ 125.

Lemma sgn_cases :
  forall x, (0 < x /\ Z.sgn x = 1) \/ (x = 0 /\ Z.sgn x = 0) \/ (x < 0 /\ Z.sgn x = -1).
Proof. intros x; destruct (Z.sgn_spec x) as [H|[H|H]]; lia. Qed.

Lemma sign_diff :
  forall P Q,
  (if negb (Z.sgn P =? Z.sgn Q) then (if Z.sgn Q <? Z.sgn P then 1 else -1)
   else Z.sgn P * Z.sgn (Z.abs P - Z.abs Q)) = Z.sgn (P - Q).
Proof.
  intros P Q.
  destruct (sgn_cases P) as [[HP EP]|[[HP EP]|[HP EP]]];
    destruct (sgn_cases Q) as [[HQ EQ]|[[HQ EQ]|[HQ EQ]]]; rewrite EP, EQ;
    cbv beta iota zeta delta [negb Z.eqb Z.ltb Z.compare Pos.compare Pos.compare_cont Pos.eqb
        CompOpp];
    try (symmetry; first [apply Z.sgn_pos | apply Z.sgn_neg]; lia).
  - rewrite (Z.abs_eq P), (Z.abs_eq Q) by lia; ring.
  - subst P Q; reflexivity.
  - rewrite (Z.abs_neq P), (Z.abs_neq Q) by lia.
    replace (- P - - Q) with (- (P - Q)) by ring; rewrite Z.sgn_opp; ring.
Qed.

Lemma sign_product_ok :
  forall f s, pword f -> pword s -> Z.abs (tc2 f) < W * W -> Z.abs (tc2 s) < W * W ->
  sign_product (times_wide (absolute f) (absolute s)) f s = Z.sgn (tc2 f * tc2 s) /\
  wide4 (times_wide (absolute f) (absolute s)) /\
  value (times_wide (absolute f) (absolute s)) = Z.abs (tc2 f * tc2 s).
Proof.
  intros f s Hf Hs Af As.
  destruct (absolute_ok f Hf) as [Pf Vf]; destruct (absolute_ok s Hs) as [Ps Vs].
  destruct (times_wide_ok _ _ Pf Ps) as [Hw V].
  rewrite Vf, Vs, <- Z.abs_mul in V.
  repeat split; try apply Hw; try exact V.
  unfold sign_product.
  destruct (forallb (fun x => x =? 0) _) eqn:Z0.
  - apply forallb_zeros, zeros_value in Z0; rewrite Z0 in V.
    symmetry in V; apply (proj1 (Z.abs_0_iff _)) in V; rewrite V; reflexivity.
  - assert (N : tc2 f * tc2 s <> 0).
    { intros E; rewrite E in V; change (Z.abs 0) with 0 in V.
      apply value_zero, forallb_zeros in V; [congruence | apply Hw]. }
    rewrite !below_tc2 by assumption.
    destruct (Z.ltb_spec (tc2 f) 0); destruct (Z.ltb_spec (tc2 s) 0); cbn [xorb]; symmetry;
      first [apply Z.sgn_pos | apply Z.sgn_neg]; nia.
Qed.

Theorem orient_grid_ok :
  forall ax ay bx by' cx cy,
  unit125 ax -> unit125 ay -> unit125 bx -> unit125 by' -> unit125 cx -> unit125 cy ->
  orient_grid ax ay bx by' cx cy =
  Z.sgn ((tc2 ax - tc2 cx) * (tc2 by' - tc2 cy) - (tc2 ay - tc2 cy) * (tc2 bx - tc2 cx)).
Proof.
  intros ax ay bx by' cx cy [Pax Aax] [Pay Aay] [Pbx Abx] [Pby Aby] [Pcx Acx] [Pcy Acy].
  unfold orient_grid; cbv zeta.
  assert (Half : HalfW * W = 2 ^ 127) by (rewrite W_pow; reflexivity).
  assert (Two : W * W = 2 ^ 128) by (rewrite W_pow; reflexivity).
  destruct (minus2_tc ax cx Pax Pcx ltac:(lia)) as [Pfx Tfx].
  destruct (minus2_tc ay cy Pay Pcy ltac:(lia)) as [Pfy Tfy].
  destruct (minus2_tc bx cx Pbx Pcx ltac:(lia)) as [Psx Tsx].
  destruct (minus2_tc by' cy Pby Pcy ltac:(lia)) as [Psy Tsy].
  destruct (sign_product_ok _ _ Pfx Psy ltac:(lia) ltac:(lia)) as [Sl [Wl Vl]].
  destruct (sign_product_ok _ _ Pfy Psx ltac:(lia) ltac:(lia)) as [Sr [Wr Vr]].
  rewrite Sl, Sr.
  rewrite compare_down_ok; [| apply words_rev, Wl | apply words_rev, Wr
                           | rewrite !length_rev, (proj2 Wl), (proj2 Wr); reflexivity].
  rewrite <- !value_top, Vl, Vr, sign_diff, Tfx, Tfy, Tsx, Tsy; reflexivity.
Qed.

(* Units: a double's significand, read, shifted onto the grid, then negated if the double is. *)
Definition units (sig exp lowest : Z) (negative : bool) : pair :=
  if sig =? 0 then Pair 0 0 else
  let shift := exp - lowest in
  negate (if shift <? 0 then Pair 0 (sig / 2 ^ (- shift))
          else if shift <? 64 then Pair ((sig / 2) / 2 ^ (63 - shift)) ((sig * 2 ^ shift) mod W)
          else Pair ((sig * 2 ^ (shift - 64)) mod W) 0) negative.

Definition scaled sig shift := if shift <? 0 then sig / 2 ^ (- shift) else sig * 2 ^ shift.

Lemma negate_tc :
  forall p n, pword p -> pvalue p < HalfW * W ->
  pword (negate p n) /\ tc2 (negate p n) = if n then - pvalue p else pvalue p.
Proof.
  intros p n Hp B; destruct (negate_ok p n Hp) as [Pn Vn]; split; [assumption|].
  pose proof (pvalue_bounds p Hp).
  apply tc2_of; [assumption | |].
  - rewrite Vn; destruct n; [rewrite <- (Z.mod_add (- pvalue p) 1 W2) by (pose proof W2_pos; lia);
      f_equal; ring | symmetry; apply Z.mod_small; lia].
  - destruct n; lia.
Qed.

Theorem units_ok :
  forall sig exp lowest n, 0 <= sig < 2 ^ 53 -> exp - lowest <= 72 ->
  pword (units sig exp lowest n) /\
  tc2 (units sig exp lowest n) = (if n then - scaled sig (exp - lowest) else scaled sig (exp -
      lowest))
  /\ 0 <= scaled sig (exp - lowest) < 2 ^ 125.
Proof.
  intros sig exp lowest n S B; unfold units, scaled.
  assert (Wp : W = 2 ^ 64) by apply W_pow.
  assert (Half : HalfW * W = 2 ^ 127) by (rewrite Wp; reflexivity).
  destruct (Z.eqb_spec sig 0) as [Z0|Z0].
  - subst sig; cbv beta iota.
    assert (S0 : (if exp - lowest <? 0 then 0 / 2 ^ (- (exp - lowest)) else 0 * 2 ^ (exp - lowest))
                 = 0) by (destruct (exp - lowest <? 0); [apply Zdiv_0_l | apply Z.mul_0_l]).
    rewrite S0.
    assert (T0 : tc2 (Pair 0 0) = 0).
    { unfold tc2, below, top_bit, pvalue; cbn [phigh plow]; rewrite Zdiv_0_l; cbn [Z.eqb negb];
        ring. }
    rewrite T0; repeat split; try (unfold word; pose proof W_gt; lia); destruct n; reflexivity.
  - set (shift := exp - lowest) in *.
    destruct (Z.ltb_spec shift 0).
    + assert (0 <= sig / 2 ^ (- shift) <= sig).
      { split; [apply Z.div_pos; [lia | apply Z.pow_pos_nonneg; lia]|].
        apply Z.div_le_upper_bound; [apply Z.pow_pos_nonneg; lia|].
        pose proof (Z.pow_pos_nonneg 2 (- shift) ltac:(lia) ltac:(lia)); nia. }
      assert (P : pword (Pair 0 (sig / 2 ^ (- shift))))
        by (unfold pword, word; cbn [phigh plow]; rewrite Wp; split; lia).
      destruct (negate_tc _ n P ltac:(unfold pvalue; cbn [phigh plow]; lia)) as [Pn Tn].
      unfold pvalue in Tn; cbn [phigh plow] in Tn.
      repeat split; try apply Pn; try lia.
      rewrite Tn; destruct n; ring.
    + assert (0 <= sig * 2 ^ shift < 2 ^ 125).
      { split; [apply Z.mul_nonneg_nonneg; [lia | apply Z.pow_nonneg; lia]|].
        replace (2 ^ 125) with (2 ^ 53 * 2 ^ 72) by reflexivity.
        pose proof (Z.pow_le_mono_r 2 shift 72 ltac:(lia) ltac:(lia)).
        pose proof (Z.pow_pos_nonneg 2 shift ltac:(lia) ltac:(lia)); nia. }
      destruct (Z.ltb_spec shift 64).
      * set (low := (sig * 2 ^ shift) mod W).
        set (high := sig / 2 / 2 ^ (63 - shift)).
        assert (Hh : high = sig * 2 ^ shift / W).
        { unfold high; rewrite Z.div_div by (lia || (apply Z.pow_pos_nonneg; lia)).
          replace (2 * 2 ^ (63 - shift)) with (2 ^ (64 - shift))
            by (rewrite <- Z.pow_succ_r by lia; f_equal; lia).
          rewrite Wp; replace (2 ^ 64) with (2 ^ (64 - shift) * 2 ^ shift)
            by (rewrite <- Z.pow_add_r by lia; f_equal; lia).
          rewrite Z.div_mul_cancel_r; [reflexivity | apply Z.pow_nonzero; lia
                                     | apply Z.pow_nonzero; lia]. }
        pose proof (Z.div_mod (sig * 2 ^ shift) W ltac:(pose proof W_gt; lia)).
        pose proof (Z.mod_pos_bound (sig * 2 ^ shift) W ltac:(pose proof W_gt; lia)).
        assert (0 <= sig * 2 ^ shift / W < W)
          by (split; [apply Z.div_pos | apply Z.div_lt_upper_bound]; rewrite ?Wp in *; lia).
        assert (P : pword (Pair high low)) by (unfold pword, word; cbn [phigh plow]; lia).
        assert (V : pvalue (Pair high low) = sig * 2 ^ shift) by (unfold pvalue; cbn; lia).
        destruct (negate_tc _ n P ltac:(lia)) as [Pn Tn].
        repeat split; try apply Pn; try lia.
        rewrite Tn, V; destruct n; ring.
      * set (high := (sig * 2 ^ (shift - 64)) mod W).
        assert (E : sig * 2 ^ shift = W * (sig * 2 ^ (shift - 64))).
        { rewrite Wp; replace shift with (64 + (shift - 64)) at 1 by lia.
          rewrite Z.pow_add_r by lia; ring. }
        assert (0 <= sig * 2 ^ (shift - 64) < W).
        { split; [apply Z.mul_nonneg_nonneg; [lia | apply Z.pow_nonneg; lia]|]. rewrite Wp in *;
            nia. }
        assert (Hh : high = sig * 2 ^ (shift - 64)) by (apply Z.mod_small; lia).
        assert (P : pword (Pair high 0))
          by (unfold pword, word; cbn [phigh plow]; split; [lia | pose proof W_gt; lia]).
        assert (V : pvalue (Pair high 0) = sig * 2 ^ shift) by (unfold pvalue; cbn; lia).
        destruct (negate_tc _ n P ltac:(lia)) as [Pn Tn].
        repeat split; try apply Pn; try lia.
        rewrite Tn, V; destruct n; ring.
Qed.

(* ---------------------------------------------------------------------------------------------- *)
(* Reading a double: grid.cpp's and whole.cpp's Read, against Flocq's decoding of the same bits.  *)

(* (bits >> 52) & 0x7FF, bits & ((1 << 52) - 1) and (bits >> 63) != 0. *)
Definition read_biased x := Z.land (Z.shiftr x 52) 2047.
Definition read_fraction x := Z.land x 4503599627370495.
Definition read_negative x := negb (Z.shiftr x 63 =? 0).

(* whole.cpp's Read: a subnormal's significand is its fraction, its place -1074; grid.cpp's, for
   normal doubles only, is the second branch. *)
Definition place_sig x := if read_biased x =? 0 then read_fraction x else read_fraction x + 2 ^ 52.
Definition place_exp x := if read_biased x =? 0 then -1074 else read_biased x - 1075.

Lemma read_biased_eq : forall x, 0 <= x -> read_biased x = (x / 2 ^ 52) mod 2 ^ 11.
Proof.
  intros x Hx; unfold read_biased; rewrite Z.shiftr_div_pow2 by lia.
  change 2047 with (Z.ones 11); rewrite Z.land_ones by lia; reflexivity.
Qed.

Lemma read_fraction_eq : forall x, 0 <= x -> read_fraction x = x mod 2 ^ 52.
Proof.
  intros x Hx; unfold read_fraction.
  change 4503599627370495 with (Z.ones 52); rewrite Z.land_ones by lia; reflexivity.
Qed.

Lemma read_negative_eq : forall x, 0 <= x < 2 ^ 64 -> read_negative x = (2 ^ 63 <=? x).
Proof.
  intros x Hx; unfold read_negative; rewrite Z.shiftr_div_pow2 by lia.
  destruct (Z.leb_spec (2 ^ 63) x).
  - assert (x / 2 ^ 63 = 1) as -> by (symmetry; apply Z.div_unique with (r := x - 2 ^ 63); lia).
    reflexivity.
  - rewrite Z.div_small by lia; reflexivity.
Qed.

Theorem read_ok :
  forall x, 0 <= x < 2 ^ 64 -> read_biased x <> 2047 ->
  B2R 53 1024 (b64_of_bits x)
  = (IZR ((if read_negative x then -1 else 1) * place_sig x) * bpow radix2 (place_exp x))%R.
Proof.
  intros x Hx Hb.
  unfold place_sig, place_exp; rewrite read_negative_eq by assumption.
  rewrite read_biased_eq, read_fraction_eq in * by lia.
  unfold b64_of_bits, binary_float_of_bits.
  change (FF2B (52 + 1) (2 ^ (11 - 1))) with (FF2B 53 1024).
  rewrite B2R_FF2B; unfold binary_float_of_bits_aux, split_bits.
  change (SpecFloat.emin (52 + 1) (2 ^ (11 - 1))) with (-1074).
  change (2 ^ 11 - 1) with 2047.
  replace (2 ^ 52 * 2 ^ 11) with (2 ^ 63) by reflexivity.
  pose proof (Z.mod_pos_bound x (2 ^ 52) ltac:(lia)).
  destruct (Z.eqb_spec ((x / 2 ^ 52) mod 2 ^ 11) 0) as [E0|E0].
  - destruct (x mod 2 ^ 52) as [|m|m] eqn:M; [| |lia]; cbn [FF2R].
    + destruct (2 ^ 63 <=? x); cbn; ring.
    + unfold F2R; cbn [Fnum Fexp].
      destruct (2 ^ 63 <=? x); cbn [cond_Zopp]; do 2 f_equal; ring.
  - destruct (Z.eqb_spec ((x / 2 ^ 52) mod 2 ^ 11) 2047) as [E1|E1]; [contradiction|].
    destruct (x mod 2 ^ 52 + 2 ^ 52) as [|m|m] eqn:M; [lia | | lia]; cbn [FF2R].
    unfold F2R; cbn [Fnum Fexp].
    replace ((x / 2 ^ 52) mod 2 ^ 11 + -1074 - 1) with ((x / 2 ^ 52) mod 2 ^ 11 - 1075) by ring.
    destruct (2 ^ 63 <=? x); cbn [cond_Zopp]; do 2 f_equal; ring.
Qed.

(* ---------------------------------------------------------------------------------------------- *)
(* whole.cpp's Units: a double's significand shifted onto the grid, in as many words as it takes. *)

Lemma base_pow : forall n, base n = 2 ^ (64 * Z.of_nat n).
Proof.
  induction n as [|n IH]; [reflexivity|].
  cbn [base]; rewrite IH, W_pow, Nat2Z.inj_succ, <- Z.pow_add_r by lia; f_equal; lia.
Qed.

(* The words the C++ sets, lowest first; the rest of its array is zero, which Trim drops. *)
Definition units_words (sig shift : Z) : list Z :=
  if shift <? 0 then [sig / 2 ^ (- shift)]
  else repeat 0 (Z.to_nat (shift / 64)) ++
       [(sig * 2 ^ (shift mod 64)) mod W;
        if shift mod 64 =? 0 then 0 else sig / 2 ^ (64 - shift mod 64)].

Definition units_whole (sig exp lowest : Z) (negative : bool) : whole :=
  if sig =? 0 then Whole [] negative
  else Whole (trim (units_words sig (exp - lowest))) negative.

Theorem units_whole_ok :
  forall sig exp lowest n, 0 <= sig < 2 ^ 53 ->
  good (units_whole sig exp lowest n) /\
  signed (units_whole sig exp lowest n)
  = (if n then - scaled sig (exp - lowest) else scaled sig (exp - lowest)).
Proof.
  intros sig exp lowest n S; unfold units_whole, scaled.
  assert (Wp : W = 2 ^ 64) by apply W_pow.
  destruct (Z.eqb_spec sig 0) as [Z0|Z0].
  - subst sig; change (0 =? 0) with true; cbv beta iota.
    unfold good, trimmed, signed; cbn [mag neg value rev hd app].
    split; [split; [apply Forall_nil | discriminate]|].
    destruct (exp - lowest <? 0); [rewrite Zdiv_0_l | rewrite Z.mul_0_l]; destruct n; reflexivity.
  - set (shift := exp - lowest).
    assert (Hw : words (units_words sig shift) /\ value (units_words sig shift)
                 = (if shift <? 0 then sig / 2 ^ (- shift) else sig * 2 ^ shift)).
    { unfold units_words; destruct (Z.ltb_spec shift 0).
      - assert (0 <= sig / 2 ^ (- shift) <= sig).
        { split; [apply Z.div_pos; [lia | apply Z.pow_pos_nonneg; lia]|].
          apply Z.div_le_upper_bound; [apply Z.pow_pos_nonneg; lia|].
          pose proof (Z.pow_pos_nonneg 2 (- shift) ltac:(lia) ltac:(lia)); nia. }
        split; [constructor; [unfold word; rewrite Wp; lia | constructor] | cbn [value]; ring].
      - set (k := shift / 64); set (b := shift mod 64).
        pose proof (Z.div_mod shift 64 ltac:(lia));
            pose proof (Z.mod_pos_bound shift 64 ltac:(lia)).
        assert (0 <= k) by (apply Z.div_pos; lia).
        assert (Lo : 0 <= (sig * 2 ^ b) mod W < W) by (apply Z.mod_pos_bound; pose proof W_gt; lia).
        assert (Hi : sig * 2 ^ b = W * (if b =? 0 then 0 else sig / 2 ^ (64 - b))
                                   + (sig * 2 ^ b) mod W
                     /\ 0 <= (if b =? 0 then 0 else sig / 2 ^ (64 - b)) < W).
        { destruct (Z.eqb_spec b 0) as [B0|B0].
          - rewrite B0, Z.pow_0_r, Z.mul_1_r, Z.mod_small by (rewrite Wp; lia). split; lia.
          - assert (D : sig * 2 ^ b / W = sig / 2 ^ (64 - b)).
            { rewrite Wp; replace (2 ^ 64) with (2 ^ (64 - b) * 2 ^ b)
                by (rewrite <- Z.pow_add_r by lia; f_equal; lia).
              apply Z.div_mul_cancel_r; apply Z.pow_nonzero; lia. }
            rewrite <- D; split; [apply Z.div_mod; pose proof W_gt; lia|].
            split; [apply Z.div_pos; [apply Z.mul_nonneg_nonneg; [lia | apply Z.pow_nonneg; lia]
                                     | pose proof W_gt; lia]|].
            apply Z.div_lt_upper_bound; [pose proof W_gt; lia|].
            rewrite Wp; replace (2 ^ 64 * 2 ^ 64) with (2 ^ 53 * 2 ^ 75) by reflexivity.
            pose proof (Z.pow_le_mono_r 2 b 75 ltac:(lia) ltac:(lia)).
            pose proof (Z.pow_pos_nonneg 2 b ltac:(lia) ltac:(lia)); nia. }
        destruct Hi as [Hi Bh].
        split.
        + apply words_app.
          * apply Forall_forall; intros w Hw; apply repeat_spec in Hw; subst; apply word_zero.
          * constructor; [exact Lo | constructor; [exact Bh | constructor]].
        + rewrite value_app, repeat_length, zeros_value by apply zeros_repeat.
          cbn [value]; rewrite base_pow, Z2Nat.id by lia.
          assert (Es : sig * 2 ^ shift = 2 ^ (64 * k) * (sig * 2 ^ b)).
          { replace shift with (64 * k + b) by lia; rewrite Z.pow_add_r by lia.
            rewrite Z.mul_assoc, (Z.mul_comm sig), <- Z.mul_assoc; reflexivity. }
          rewrite Es.
          set (lo := (sig * 2 ^ b) mod W) in *.
          set (hi := if b =? 0 then 0 else sig / 2 ^ (64 - b)) in *.
          rewrite Hi; clearbody lo hi k b; ring. }
    destruct Hw as [Hw V]; destruct (trim_ok _ Hw) as [Hw' [T [V' _]]].
    split; [split; assumption|].
    unfold signed; cbn [mag neg]; rewrite V', V; destruct n; reflexivity.
Qed.

(* ---------------------------------------------------------------------------------------------- *)
(* The signs, against the real determinants of the doubles given.                                 *)

Definition coordinate x := B2R 53 1024 (b64_of_bits x).

(* A double the code reads: its bits finite. *)
Definition finite_bits x := 0 <= x < 2 ^ 64 /\ read_biased x <> 2047.

(* The grid the code picks lies at or below every set bit: Lowest takes the least place plus
   std::countr_zero, which 2 to its power divides the significand. *)
Definition under x L :=
  place_sig x = 0 \/ exists t, 0 <= t /\ (2 ^ t | place_sig x) /\ L <= place_exp x + t.

Definition signum (r : R) : Z :=
  match Rcompare r 0 with Lt => -1 | Eq => 0 | Gt => 1 end.

Lemma place_sig_bound : forall x, finite_bits x -> 0 <= place_sig x < 2 ^ 53.
Proof.
  intros x [Hx _]; unfold place_sig; rewrite read_fraction_eq by lia.
  pose proof (Z.mod_pos_bound x (2 ^ 52) ltac:(lia)).
  destruct (read_biased x =? 0); lia.
Qed.

Lemma IZR_pow2 : forall e, 0 <= e -> IZR (2 ^ e) = bpow radix2 e.
Proof. intros [|e|e] He; [reflexivity | reflexivity | lia]. Qed.

Lemma scaled_place :
  forall sig exp L, 0 <= sig -> (sig = 0 \/ exists t, 0 <= t /\ (2 ^ t | sig) /\ L <= exp + t) ->
  (IZR (scaled sig (exp - L)) * bpow radix2 L = IZR sig * bpow radix2 exp)%R.
Proof.
  intros sig exp L S U; unfold scaled.
  destruct (Z.ltb_spec (exp - L) 0).
  - destruct U as [Z0|[t [Ht [[q Hq] Hl]]]].
    + subst sig; rewrite Zdiv_0_l; cbn; ring.
    + assert (E : sig = (q * 2 ^ (t - (L - exp))) * 2 ^ (L - exp)).
      { rewrite Hq, <- Z.mul_assoc, <- Z.pow_add_r by lia; f_equal; f_equal; lia. }
      replace (- (exp - L)) with (L - exp) by ring.
      rewrite E, Z.div_mul by (apply Z.pow_nonzero; lia).
      rewrite !mult_IZR, (IZR_pow2 (L - exp)) by lia.
      rewrite (Rmult_assoc _ (bpow radix2 (L - exp)) (bpow radix2 exp)), <- bpow_plus.
      replace (L - exp + exp) with L by ring; reflexivity.
  - rewrite mult_IZR, (IZR_pow2 (exp - L)) by lia.
    rewrite (Rmult_assoc _ (bpow radix2 (exp - L)) (bpow radix2 L)), <- bpow_plus.
    replace (exp - L + L) with exp by ring; reflexivity.
Qed.

Lemma coordinate_scaled :
  forall x L, finite_bits x -> under x L ->
  coordinate x =
  (IZR (if read_negative x then - scaled (place_sig x) (place_exp x - L)
        else scaled (place_sig x) (place_exp x - L)) * bpow radix2 L)%R.
Proof.
  intros x L F U; unfold coordinate; rewrite read_ok by apply F.
  pose proof (place_sig_bound x F).
  pose proof (scaled_place (place_sig x) (place_exp x) L ltac:(lia) U) as E.
  destruct (read_negative x); rewrite mult_IZR.
  - rewrite opp_IZR; change (IZR (-1)) with (- 1)%R.
    replace (- IZR (scaled (place_sig x) (place_exp x - L)) * bpow radix2 L)%R
      with (- (IZR (scaled (place_sig x) (place_exp x - L)) * bpow radix2 L))%R by ring.
    rewrite E; ring.
  - change (IZR 1) with 1%R; rewrite E; ring.
Qed.

Lemma signum_scale : forall z k, (0 < k)%R -> signum (IZR z * k) = Z.sgn z.
Proof.
  intros z k K; unfold signum.
  destruct (Z.compare_spec z 0) as [E|E|E].
  - subst; rewrite Rmult_0_l, Rcompare_Eq; reflexivity.
  - rewrite Rcompare_Lt; [rewrite Z.sgn_neg by lia; reflexivity|].
    replace 0%R with (0 * k)%R by ring; apply Rmult_lt_compat_r; [exact K|].
    apply IZR_lt; exact E.
  - rewrite Rcompare_Gt; [rewrite Z.sgn_pos by lia; reflexivity|].
    replace 0%R with (0 * k)%R by ring; apply Rmult_lt_compat_r; [exact K|].
    apply IZR_lt; exact E.
Qed.

Definition orient_real (ax ay bx by' cx cy : R) :=
  ((ax - cx) * (by' - cy) - (ay - cy) * (bx - cx))%R.

Definition incircle_real (ax ay bx by' cx cy : R) :=
  ((ax * ax + ay * ay) * (bx * cy - by' * cx) + (bx * bx + by' * by') * (cx * ay - cy * ax)
   + (cx * cx + cy * cy) * (ax * by' - ay * bx))%R.

Lemma orient_scaled :
  forall a b c d e f k,
  orient_real (IZR a * k) (IZR b * k) (IZR c * k) (IZR d * k) (IZR e * k) (IZR f * k)
  = (IZR ((a - e) * (d - f) - (b - f) * (c - e)) * (k * k))%R.
Proof. intros; unfold orient_real; rewrite !minus_IZR, !mult_IZR, !minus_IZR; ring. Qed.

Lemma incircle_scaled :
  forall a b c d e f g h k,
  incircle_real (IZR a * k - IZR g * k) (IZR b * k - IZR h * k) (IZR c * k - IZR g * k)
                (IZR d * k - IZR h * k) (IZR e * k - IZR g * k) (IZR f * k - IZR h * k)
  = (IZR (incircle (a - g) (b - h) (c - g) (d - h) (e - g) (f - h)) * (k * k * (k * k)))%R.
Proof.
  intros; unfold incircle_real, incircle.
  rewrite !plus_IZR, !mult_IZR, !minus_IZR, !plus_IZR, !mult_IZR, !minus_IZR; ring.
Qed.

Lemma scaled_small :
  forall sig shift k, 0 <= sig < 2 ^ 53 -> shift <= k -> 0 <= k ->
  0 <= scaled sig shift < 2 ^ (53 + k).
Proof.
  intros sig shift k S B K; unfold scaled; destruct (Z.ltb_spec shift 0).
  - assert (2 ^ 53 <= 2 ^ (53 + k)) by (apply Z.pow_le_mono_r; lia).
    split; [apply Z.div_pos; [lia | apply Z.pow_pos_nonneg; lia]|].
    assert (sig / 2 ^ (- shift) <= sig).
    { apply Z.div_le_upper_bound; [apply Z.pow_pos_nonneg; lia|].
      pose proof (Z.pow_pos_nonneg 2 (- shift) ltac:(lia) ltac:(lia)); nia. }
    lia.
  - rewrite Z.pow_add_r by lia.
    pose proof (Z.pow_le_mono_r 2 shift k ltac:(lia) ltac:(lia)).
    pose proof (Z.pow_pos_nonneg 2 shift ltac:(lia) ltac:(lia)).
    split; [apply Z.mul_nonneg_nonneg; lia | nia].
Qed.

(* A coordinate read onto the grid: its units, two words of them, are the double over 2^L. *)
Lemma read_unit :
  forall x L k, finite_bits x -> under x L -> place_exp x - L <= k -> 0 <= k <= 72 ->
  let u := units (place_sig x) (place_exp x) L (read_negative x) in
  pword u /\ Z.abs (tc2 u) < 2 ^ (53 + k) /\ coordinate x = (IZR (tc2 u) * bpow radix2 L)%R.
Proof.
  intros x L k F U B K u.
  pose proof (place_sig_bound x F) as S.
  destruct (units_ok (place_sig x) (place_exp x) L (read_negative x) S ltac:(lia)) as [P [T _]].
  pose proof (scaled_small (place_sig x) (place_exp x - L) k S B ltac:(lia)).
  refine (conj P (conj _ _)).
  - fold u in T; rewrite T; destruct (read_negative x); lia.
  - rewrite (coordinate_scaled x L) by assumption; fold u in T; rewrite T; reflexivity.
Qed.

Lemma bpow_square : forall L, (0 < bpow radix2 L * bpow radix2 L)%R.
Proof. intros L; apply Rmult_lt_0_compat; apply bpow_gt_0. Qed.

(* OnGrid's Count: a double's significand shifted onto the grid in one word, as a std::int64_t:
   sig << shift or sig >> -shift, then negated if the double is; exact below 2^61. *)
Definition count (sig exp lowest : Z) (negative : bool) : Z :=
  if sig =? 0 then 0 else
  let shift := exp - lowest in
  let whole := if 0 <=? shift then (sig * 2 ^ shift) mod W else sig / 2 ^ (- shift) in
  if negative then - whole else whole.

Lemma read_count :
  forall x L, finite_bits x -> under x L -> place_exp x - L <= 8 ->
  let c := count (place_sig x) (place_exp x) L (read_negative x) in
  Z.abs c < 2 ^ 61 /\ coordinate x = (IZR c * bpow radix2 L)%R.
Proof.
  intros x L F U B c.
  pose proof (place_sig_bound x F) as S.
  pose proof (scaled_small (place_sig x) (place_exp x - L) 8 S B ltac:(lia)) as Sc.
  assert (E : c = if read_negative x then - scaled (place_sig x) (place_exp x - L)
                  else scaled (place_sig x) (place_exp x - L)).
  { unfold c, count, scaled; cbv zeta.
    destruct (Z.eqb_spec (place_sig x) 0) as [Z0|Z0].
    - rewrite Z0, Zdiv_0_l, Z.mul_0_l; destruct (place_exp x - L <? 0), (read_negative x);
        reflexivity.
    - destruct (Z.leb_spec 0 (place_exp x - L)); destruct (Z.ltb_spec (place_exp x - L) 0);
        try lia.
      rewrite Z.mod_small; [reflexivity|].
      split; [apply Z.mul_nonneg_nonneg; [lia | apply Z.pow_nonneg; lia]|].
      unfold scaled in Sc; destruct (Z.ltb_spec (place_exp x - L) 0); [lia|].
      rewrite W_pow; lia. }
  split.
  - rewrite E; destruct (read_negative x); lia.
  - rewrite (coordinate_scaled x L) by assumption; rewrite E; reflexivity.
Qed.

(* OrientGrid: where every coordinate's top bit lies at most 124 places above the grid, as its
   check of the span asks, its sign is the sign of the determinant of the doubles. *)
Theorem orient_grid_correct :
  forall ax ay bx by' cx cy L,
  finite_bits ax -> finite_bits ay -> finite_bits bx -> finite_bits by' -> finite_bits cx ->
  finite_bits cy ->
  under ax L -> under ay L -> under bx L -> under by' L -> under cx L -> under cy L ->
  place_exp ax - L <= 72 -> place_exp ay - L <= 72 -> place_exp bx - L <= 72 ->
  place_exp by' - L <= 72 -> place_exp cx - L <= 72 -> place_exp cy - L <= 72 ->
  let u x := units (place_sig x) (place_exp x) L (read_negative x) in
  orient_grid (u ax) (u ay) (u bx) (u by') (u cx) (u cy) =
  signum (orient_real (coordinate ax) (coordinate ay) (coordinate bx) (coordinate by')
                      (coordinate cx) (coordinate cy)).
Proof.
  intros ax ay bx by' cx cy L Fa Fb Fc Fd Fe Ff Ua Ub Uc Ud Ue Uf Ba Bb Bc Bd Be Bf u.
  destruct (read_unit ax L 72 Fa Ua Ba ltac:(lia)) as [Pa [Aa Ca]].
  destruct (read_unit ay L 72 Fb Ub Bb ltac:(lia)) as [Pb [Ab Cb]].
  destruct (read_unit bx L 72 Fc Uc Bc ltac:(lia)) as [Pc [Ac Cc]].
  destruct (read_unit by' L 72 Fd Ud Bd ltac:(lia)) as [Pd [Ad Cd]].
  destruct (read_unit cx L 72 Fe Ue Be ltac:(lia)) as [Pe [Ae Ce]].
  destruct (read_unit cy L 72 Ff Uf Bf ltac:(lia)) as [Pf [Af Cf]].
  rewrite orient_grid_ok by (split; assumption).
  rewrite Ca, Cb, Cc, Cd, Ce, Cf, orient_scaled, signum_scale by apply bpow_square.
  reflexivity.
Qed.

(* OnGrid: where every coordinate's top bit lies at most 60 places above the grid. *)
Theorem on_grid_correct :
  forall ax ay bx by' cx cy dx dy L,
  finite_bits ax -> finite_bits ay -> finite_bits bx -> finite_bits by' -> finite_bits cx ->
  finite_bits cy -> finite_bits dx -> finite_bits dy ->
  under ax L -> under ay L -> under bx L -> under by' L -> under cx L -> under cy L ->
  under dx L -> under dy L ->
  place_exp ax - L <= 8 -> place_exp ay - L <= 8 -> place_exp bx - L <= 8 ->
  place_exp by' - L <= 8 -> place_exp cx - L <= 8 -> place_exp cy - L <= 8 ->
  place_exp dx - L <= 8 -> place_exp dy - L <= 8 ->
  let u x := count (place_sig x) (place_exp x) L (read_negative x) in
  on_grid (u ax) (u ay) (u bx) (u by') (u cx) (u cy) (u dx) (u dy) =
  signum (incircle_real (coordinate ax - coordinate dx) (coordinate ay - coordinate dy)
                        (coordinate bx - coordinate dx) (coordinate by' - coordinate dy)
                        (coordinate cx - coordinate dx) (coordinate cy - coordinate dy)).
Proof.
  intros ax ay bx by' cx cy dx dy L Fa Fb Fc Fd Fe Ff Fg Fh Ua Ub Uc Ud Ue Uf Ug Uh
    Ba Bb Bc Bd Be Bf Bg Bh u.
  destruct (read_count ax L Fa Ua Ba) as [Aa Ca].
  destruct (read_count ay L Fb Ub Bb) as [Ab Cb].
  destruct (read_count bx L Fc Uc Bc) as [Ac Cc].
  destruct (read_count by' L Fd Ud Bd) as [Ad Cd].
  destruct (read_count cx L Fe Ue Be) as [Ae Ce].
  destruct (read_count cy L Ff Uf Bf) as [Af Cf].
  destruct (read_count dx L Fg Ug Bg) as [Ag Cg].
  destruct (read_count dy L Fh Uh Bh) as [Ah Ch].
  rewrite on_grid_ok by assumption.
  rewrite Ca, Cb, Cc, Cd, Ce, Cf, Cg, Ch, incircle_scaled, signum_scale
    by (apply Rmult_lt_0_compat; apply bpow_square).
  reflexivity.
Qed.

Lemma read_whole :
  forall x L, finite_bits x -> under x L ->
  let w := units_whole (place_sig x) (place_exp x) L (read_negative x) in
  good w /\ coordinate x = (IZR (signed w) * bpow radix2 L)%R.
Proof.
  intros x L F U w.
  pose proof (place_sig_bound x F) as S.
  destruct (units_whole_ok (place_sig x) (place_exp x) L (read_negative x) S) as [G T].
  split; [exact G|]; rewrite (coordinate_scaled x L) by assumption; fold w in T; rewrite T;
  reflexivity.
Qed.

(* OrientWhole, for every finite double. *)
Theorem orient_whole_correct :
  forall ax ay bx by' cx cy L,
  finite_bits ax -> finite_bits ay -> finite_bits bx -> finite_bits by' -> finite_bits cx ->
  finite_bits cy ->
  under ax L -> under ay L -> under bx L -> under by' L -> under cx L -> under cy L ->
  let w x := units_whole (place_sig x) (place_exp x) L (read_negative x) in
  orient_whole (w ax) (w ay) (w bx) (w by') (w cx) (w cy) =
  signum (orient_real (coordinate ax) (coordinate ay) (coordinate bx) (coordinate by')
                      (coordinate cx) (coordinate cy)).
Proof.
  intros ax ay bx by' cx cy L Fa Fb Fc Fd Fe Ff Ua Ub Uc Ud Ue Uf w.
  destruct (read_whole ax L Fa Ua) as [Ga Ca]; destruct (read_whole ay L Fb Ub) as [Gb Cb].
  destruct (read_whole bx L Fc Uc) as [Gc Cc]; destruct (read_whole by' L Fd Ud) as [Gd Cd].
  destruct (read_whole cx L Fe Ue) as [Ge Ce]; destruct (read_whole cy L Ff Uf) as [Gf Cf].
  rewrite orient_whole_ok by assumption.
  rewrite Ca, Cb, Cc, Cd, Ce, Cf, orient_scaled, signum_scale by apply bpow_square.
  reflexivity.
Qed.

(* CircleWhole, for every finite double. *)
Theorem circle_whole_correct :
  forall ax ay bx by' cx cy dx dy L,
  finite_bits ax -> finite_bits ay -> finite_bits bx -> finite_bits by' -> finite_bits cx ->
  finite_bits cy -> finite_bits dx -> finite_bits dy ->
  under ax L -> under ay L -> under bx L -> under by' L -> under cx L -> under cy L ->
  under dx L -> under dy L ->
  let w x := units_whole (place_sig x) (place_exp x) L (read_negative x) in
  circle_whole (w ax) (w ay) (w bx) (w by') (w cx) (w cy) (w dx) (w dy) =
  signum (incircle_real (coordinate ax - coordinate dx) (coordinate ay - coordinate dy)
                        (coordinate bx - coordinate dx) (coordinate by' - coordinate dy)
                        (coordinate cx - coordinate dx) (coordinate cy - coordinate dy)).
Proof.
  intros ax ay bx by' cx cy dx dy L Fa Fb Fc Fd Fe Ff Fg Fh Ua Ub Uc Ud Ue Uf Ug Uh w.
  destruct (read_whole ax L Fa Ua) as [Ga Ca]; destruct (read_whole ay L Fb Ub) as [Gb Cb].
  destruct (read_whole bx L Fc Uc) as [Gc Cc]; destruct (read_whole by' L Fd Ud) as [Gd Cd].
  destruct (read_whole cx L Fe Ue) as [Ge Ce]; destruct (read_whole cy L Ff Uf) as [Gf Cf].
  destruct (read_whole dx L Fg Ug) as [Gg Cg]; destruct (read_whole dy L Fh Uh) as [Gh Ch].
  rewrite circle_whole_ok by assumption.
  rewrite Ca, Cb, Cc, Cd, Ce, Cf, Cg, Ch, incircle_scaled, signum_scale
    by (apply Rmult_lt_0_compat; apply bpow_square).
  reflexivity.
Qed.

Print Assumptions multiply_halves.
Print Assumptions orient_grid_correct.
Print Assumptions on_grid_correct.
Print Assumptions orient_whole_correct.
Print Assumptions circle_whole_correct.
