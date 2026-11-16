# Two random walks meeting before a frozen trap: chains on at most four states

Yitzchak Shmalo, 16 November 2026.

Two independent rate-one random walks X and Y on a finite reversible Markov chain and a frozen point Z all start from the stationary law pi. Addario-Berry, Oliveira, Peres and Sousi (Electron. J. Probab. 21, 2016, Paper No. 42; arXiv:1506.07845) conjectured that X and Y meet before either of them hits Z with probability at least 1/3. The paper proves the stronger bound 1/3 + (2/3) sum_x pi(x)^3 for every reversible chain on three or four states, so the conjecture holds for chains on at most four states. Equivalently, unless the three points start at one state, X and Y meet first with probability greater than 1/3. On three states the proof is a closed formula and the refined constant is approached at every pi; on four states it is a polynomial with positive coefficients, built by two independent programs. The refined inequality is conjectured for all finite reversible chains. The general case of the original conjecture remains open.

## Files

- `main.tex`, `main.pdf`: the paper.
- `CREDITS.md`: which statements are known and which are new, and whose work they rest on.
- `code/meeting.py`: exact meeting probabilities for small chains by solving the first-step equations on ordered pairs in rational arithmetic; reproduces the example values of Remark 5.3.
- `code/n3.py`: checks the three-state identity (Proposition 3.1), the closed form of D against `meeting.py`, strictness, and the sharpness family.
- `code/n4_cramer.py`: builds the four-state certificate R of Section 4 with SymPy on the six support graphs, checks that every coefficient is positive and that N = R/Q at random integer points.
- `code/n4_integer.py`: an independent construction of R with plain integer polynomials (no SymPy); with `--compare` it checks agreement with `n4_cramer.py` coefficient by coefficient.
- `code/output.txt`: the output of all four programs.

Every program exits with status 1 if any check fails. They need Python 3 and SymPy (`n4_integer.py` needs only `meeting.py`). The K4 support takes several minutes and about 1 GB in `n4_cramer.py`, and longer in `n4_integer.py`.

## License

The paper is licensed CC BY 4.0 and the code MIT.
