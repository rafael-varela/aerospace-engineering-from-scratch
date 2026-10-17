# P01.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Row 2 is twice row 1, so one row will vanish; the pivot columns are the first two. A special solution has a 1 in its own free column and the negated rref entries in the pivot rows. For (d), the solvability condition is that b is orthogonal to the left null space.
2. For (a), look at the free rows of the special solutions: they form an identity block, so no combination of them can vanish; for spanning, subtract from any null vector the combination of special solutions that matches its free components and show what is left is zero. For (d), apply dim N(A) = n − r and dim N(Aᵀ) = m − r with n = m_members + r_reactions and m = 2j, then subtract.
3. The achievable torques are the column space of B; the impossible direction is the left null space; a null vector of B is a set of wheel torques that cancel. Dropping a wheel deletes a column: check the rank of what is left.
4. rref scales each pivot row to 1 and never asks how small the pivot was before scaling; the singular values do. Forces scale as 1/σ_min. The tolerance at which the rank drops compares the pivot of the diagonal's vertical component (≈ y₄/2.5) with the largest entry (1).
5. The mode's non-zero entries tell you which joints move; a repair must connect a moving joint to a joint that does not move *in a direction with a component along the motion*. A member between two joints that move identically stretches nothing and adds only a self-stress state.
6. Separate three questions: is the structure a first-order mechanism (rank), how far is it from one (σ_min, condition number), and is the force path acceptable (stiffness analysis, buckling). Each statement conflates two of them.
7. A vector **y**⊥ with **y**⊥ᵀ H = 0 kills the load term in **y**⊥ᵀ (H **f** + **e**) and leaves **y**⊥ᵀ **e**. With one parity relation the signatures are scalars and any two non-zero scalars are "parallel"; with two relations they are 2-vectors and can point in different directions.
