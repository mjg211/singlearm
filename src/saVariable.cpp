#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix saVariable(int J, double pi0, double pi1, double alpha,
                         double beta, int v, int Nmin, int Nmax, int equal_n) {
  NumericMatrix dbinomial_pi0(Nmax, Nmax), dbinomial_pi1(Nmax, Nmax);
  NumericMatrix pbinomial_pi0(Nmax, Nmax), pbinomial_pi1(Nmax, Nmax);
  for (int j = 1; j <= Nmax - 1; j++) {
    for (int i = 0; i <= j; i++) {
      dbinomial_pi0(i, j) = R::dbinom(i, j, pi0, 0);
      dbinomial_pi1(i, j) = R::dbinom(i, j, pi1, 0);
      pbinomial_pi0(i, j) = R::pbinom(i, j, pi0, 1, 0);
      pbinomial_pi1(i, j) = R::pbinom(i, j, pi1, 1, 0);
    }
  }
  int counter = 0;
  //int interrupt = 0;
  NumericMatrix feasible_designs(10000000, (J == 1 ? 2*v + 2 : 4*v + 4));
  if (J == 1) {
    NumericVector R_pi0(v), R_pi1(v);
    if (v == 2) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (sum(R_pi1) < 2 - beta*v) {
              break;
            }
            if (mean(R_pi0) <= alpha) {
              feasible_designs(counter, _) = NumericVector::create(n[0], n[1], a1, a2,
                               mean(R_pi0), mean(R_pi1));
              counter++;
            }
          }
        }
      }
    }
    else if (v == 3) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (sum(R_pi1) < 3 - beta*v) {
                break;
              }
              if (mean(R_pi0) <= alpha) {
                feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], a1, a2, a3,
                                 mean(R_pi0), mean(R_pi1));
                counter++;
              }
            }
          }
        }
      }
    }
    else if (v == 4) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (R_pi1[0] + R_pi1[1] + R_pi1[2] < 3 - beta*v) {
                break;
              }
              for (int a4 = a3; a4 < n[3]; a4++) {
                R_pi0[3]                 = 1 - pbinomial_pi0(a4, n[3]);
                R_pi1[3]                 = 1 - pbinomial_pi1(a4, n[3]);
                if (sum(R_pi1) < 4 - beta*v) {
                  break;
                }
                if (mean(R_pi0) <= alpha) {
                  feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], n[3], a1, a2, a3,
                                   a4, mean(R_pi0), mean(R_pi1));
                  counter++;
                }
              }
            }
          }
        }
      }
    }
    else if (v == 5) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (R_pi1[0] + R_pi1[1] + R_pi1[2] < 3 - beta*v) {
                break;
              }
              for (int a4 = a3; a4 < n[3]; a4++) {
                R_pi0[3]                 = 1 - pbinomial_pi0(a4, n[3]);
                R_pi1[3]                 = 1 - pbinomial_pi1(a4, n[3]);
                if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] < 4 - beta*v) {
                  break;
                }
                for (int a5 = a4; a5 < n[4]; a5++) {
                  R_pi0[4]                 = 1 - pbinomial_pi0(a5, n[4]);
                  R_pi1[4]                 = 1 - pbinomial_pi1(a5, n[4]);
                  if (sum(R_pi1) < 5 - beta*v) {
                    break;
                  }
                  if (mean(R_pi0) <= alpha) {
                    feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], n[3], n[4],
                                     a1, a2, a3, a4, a5, mean(R_pi0), mean(R_pi1));
                    counter++;
                  }
                }
              }
            }
          }
        }
      }
    }
    else if (v == 6) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (R_pi1[0] + R_pi1[1] + R_pi1[2] < 3 - beta*v) {
                break;
              }
              for (int a4 = a3; a4 < n[3]; a4++) {
                R_pi0[3]                 = 1 - pbinomial_pi0(a4, n[3]);
                R_pi1[3]                 = 1 - pbinomial_pi1(a4, n[3]);
                if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] < 4 - beta*v) {
                  break;
                }
                for (int a5 = a4; a5 < n[4]; a5++) {
                  R_pi0[4]                 = 1 - pbinomial_pi0(a5, n[4]);
                  R_pi1[4]                 = 1 - pbinomial_pi1(a5, n[4]);
                  if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] < 5 - beta*v) {
                    break;
                  }
                  for (int a6 = a5; a6 < n[5]; a6++) {
                    R_pi0[5]                 = 1 - pbinomial_pi0(a6, n[5]);
                    R_pi1[5]                 = 1 - pbinomial_pi1(a6, n[5]);
                    if (sum(R_pi1) < 6 - beta*v) {
                      break;
                    }
                    if (mean(R_pi0) <= alpha) {
                      feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], n[3], n[4],
                                       n[5], a1, a2, a3, a4, a5,a6,  mean(R_pi0), mean(R_pi1));
                      counter++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else if (v == 7) {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (R_pi1[0] + R_pi1[1] + R_pi1[2] < 3 - beta*v) {
                break;
              }
              for (int a4 = a3; a4 < n[3]; a4++) {
                R_pi0[3]                 = 1 - pbinomial_pi0(a4, n[3]);
                R_pi1[3]                 = 1 - pbinomial_pi1(a4, n[3]);
                if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] < 4 - beta*v) {
                  break;
                }
                for (int a5 = a4; a5 < n[4]; a5++) {
                  R_pi0[4]                 = 1 - pbinomial_pi0(a5, n[4]);
                  R_pi1[4]                 = 1 - pbinomial_pi1(a5, n[4]);
                  if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] < 5 - beta*v) {
                    break;
                  }
                  for (int a6 = a5; a6 < n[5]; a6++) {
                    R_pi0[5]                 = 1 - pbinomial_pi0(a6, n[5]);
                    R_pi1[5]                 = 1 - pbinomial_pi1(a6, n[5]);
                    if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] + R_pi1[5] < 6 - beta*v) {
                      break;
                    }
                    for (int a7 = a6; a7 < n[6]; a7++) {
                      R_pi0[6]                 = 1 - pbinomial_pi0(a7, n[6]);
                      R_pi1[6]                 = 1 - pbinomial_pi1(a7, n[6]);
                      if (sum(R_pi1) < 7 - beta*v) {
                        break;
                      }
                      if (mean(R_pi0) <= alpha) {
                        feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], n[3], n[4],
                                         n[5], n[6], a1, a2, a3, a4, a5, a6, a7, mean(R_pi0), mean(R_pi1));
                        counter++;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else {
      for (int n1 = Nmin; n1 <= Nmax - v + 1; n1++) {
        IntegerVector n = seq(n1, n1 + v - 1);
        for (int a1 = 0; a1 < n[0]; a1++) {
          R_pi0[0]                 = 1 - pbinomial_pi0(a1, n[0]);
          R_pi1[0]                 = 1 - pbinomial_pi1(a1, n[0]);
          if (R_pi1[0] < 1 - beta*v) {
            break;
          }
          for (int a2 = a1; a2 < n[1]; a2++) {
            R_pi0[1]                 = 1 - pbinomial_pi0(a2, n[1]);
            R_pi1[1]                 = 1 - pbinomial_pi1(a2, n[1]);
            if (R_pi1[0] + R_pi1[1] < 2 - beta*v) {
              break;
            }
            for (int a3 = a2; a3 < n[2]; a3++) {
              R_pi0[2]                 = 1 - pbinomial_pi0(a3, n[2]);
              R_pi1[2]                 = 1 - pbinomial_pi1(a3, n[2]);
              if (R_pi1[0] + R_pi1[1] + R_pi1[2] < 3 - beta*v) {
                break;
              }
              for (int a4 = a3; a4 < n[3]; a4++) {
                R_pi0[3]                 = 1 - pbinomial_pi0(a4, n[3]);
                R_pi1[3]                 = 1 - pbinomial_pi1(a4, n[3]);
                if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] < 4 - beta*v) {
                  break;
                }
                for (int a5 = a4; a5 < n[4]; a5++) {
                  R_pi0[4]                 = 1 - pbinomial_pi0(a5, n[4]);
                  R_pi1[4]                 = 1 - pbinomial_pi1(a5, n[4]);
                  if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] < 5 - beta*v) {
                    break;
                  }
                  for (int a6 = a5; a6 < n[5]; a6++) {
                    R_pi0[5]                 = 1 - pbinomial_pi0(a6, n[5]);
                    R_pi1[5]                 = 1 - pbinomial_pi1(a6, n[5]);
                    if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] + R_pi1[5] < 6 - beta*v) {
                      break;
                    }
                    for (int a7 = a6; a7 < n[6]; a7++) {
                      R_pi0[6]                 = 1 - pbinomial_pi0(a7, n[6]);
                      R_pi1[6]                 = 1 - pbinomial_pi1(a7, n[6]);
                      if (R_pi1[0] + R_pi1[1] + R_pi1[2] + R_pi1[3] + R_pi1[4] + R_pi1[5] + R_pi1[6] < 7 - beta*v) {
                        break;
                      }
                      for (int a8 = a7; a8 < n[7]; a8++) {
                        R_pi0[7]                 = 1 - pbinomial_pi0(a8, n[7]);
                        R_pi1[7]                 = 1 - pbinomial_pi1(a8, n[7]);
                        if (sum(R_pi1) < 8 - beta*v) {
                          break;
                        }
                        if (mean(R_pi0) <= alpha) {
                          feasible_designs(counter, _) = NumericVector::create(n[0], n[1], n[2], n[3], n[4],
                                           n[5], n[6], n[7], a1, a2, a3, a4, a5, a6, a7, a8, mean(R_pi0), mean(R_pi1));
                          counter++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else {
    NumericMatrix A1_pi0(v, v), A1_pi1(v, v);
    if (v == 2) {
      for (int n11 = 1; n11 <= Nmax - 2*v + 1; n11++) {
        IntegerVector n1 = seq(n11, n11 + v - 1);
        for (int a11 = 0; a11 < n1[0]; a11++) {
          A1_pi0(0, _) = rep(pbinomial_pi0(a11, n1[0]), v);
          A1_pi1(0, _) = rep(pbinomial_pi1(a11, n1[0]), v);
          if (sum(A1_pi1(0, _)) >= beta*pow(v, 2)) {
            break;
          }
          for (int a12 = a11; a12 < n1[1]; a12++) {
            A1_pi0(1, _) = rep(pbinomial_pi0(a12, n1[1]), v);
            A1_pi1(1, _) = rep(pbinomial_pi1(a12, n1[1]), v);
            if (sum(A1_pi1) >= beta*pow(v, 2)) {
              break;
            }
            NumericVector a1 = NumericVector::create(a11, a12);
            for (int n21 = max(1, Nmin - n11 - v + 1); n21 <= Nmax - n11 - v; n21++) {
              IntegerVector n2 = seq(n21, n21 + v - 1);
              for (int a21 = a11 + 1; a21 < n1[0] + n2[0]; a21++) {
                NumericMatrix A2_pi0(v, v), A2_pi1(v, v);
                for (int i = 0; i <= v - 1; i++) {
                  for (int s1 = a1[i] + 1; s1 <= min(a21, n1[i]); s1++) {
                    A2_pi0(i, 0) += dbinomial_pi0(s1, n1[i])*pbinomial_pi0(a21 - s1, n2[0]);
                    A2_pi1(i, 0) += dbinomial_pi1(s1, n1[i])*pbinomial_pi1(a21 - s1, n2[0]);
                  }
                }
                if (sum(A1_pi1) + sum(A2_pi1) >= beta*pow(v, 2)) {
                  break;
                }
                for (int a22 = max(a12 + 1, a21); a22 < n1[0] + n2[1]; a22++) {
                  for (int i = 0; i <= v - 1; i++) {
                    A2_pi0(i, 1) = 0;
                    A2_pi1(i, 1) = 0;
                    for (int s1 = a1[i] + 1; s1 <= min(a22, n1[i]); s1++) {
                      A2_pi0(i, 1) += dbinomial_pi0(s1, n1[i])*pbinomial_pi0(a22 - s1, n2[1]);
                      A2_pi1(i, 1) += dbinomial_pi1(s1, n1[i])*pbinomial_pi1(a22 - s1, n2[1]);
                    }
                  }
                  if (sum(A1_pi1) + sum(A2_pi1) >= beta*pow(v, 2)) {
                    break;
                  }
                  if (mean(A1_pi0) + mean(A2_pi0) > 1 - alpha) {
                    double SESS_pi0 = 0;
                    double SESS_pi1 = 0;
                    for (int i = 1; i <= v; i++) {
                      for (int j = 1; j <= v; j++) {
                        SESS_pi0 += n1[i - 1] + (1 - A1_pi0(i - 1, j - 1))*n2[j - 1];
                        SESS_pi1 += n1[i - 1] + (1 - A1_pi1(i - 1, j - 1))*n2[j - 1];
                      }
                    }
                    feasible_designs(counter, _) = NumericVector::create(n1[0], n1[1],
                                     n2[0], n2[1],
                                              a11, a12,
                                              a21, a22,
                                              1 - mean(A1_pi0) - mean(A2_pi0),
                                              1 - mean(A1_pi1) - mean(A2_pi1), SESS_pi0/pow(v, 2), SESS_pi1/pow(v, 2));
                    counter++;
                  }
                }
              }
            }
          }
        }
      }
    }
    else if (v == 3) {
      for (int n11 = 1; n11 <= Nmax - 2*v + 1; n11++) {
        IntegerVector n1 = seq(n11, n11 + v - 1);
        for (int a11 = 0; a11 < n1[0]; a11++) {
          A1_pi0(0, _) = rep(pbinomial_pi0(a11, n1[0]), v);
          A1_pi1(0, _) = rep(pbinomial_pi1(a11, n1[0]), v);
          if (sum(A1_pi1(0, _)) >= beta*pow(v, 2)) {
            break;
          }
          for (int a12 = a11; a12 < n1[1]; a12++) {
            A1_pi0(1, _) = rep(pbinomial_pi0(a12, n1[1]), v);
            A1_pi1(1, _) = rep(pbinomial_pi1(a12, n1[1]), v);
            if (sum(A1_pi1(0, _)) + sum(A1_pi1(1, _)) >= beta*pow(v, 2)) {
              break;
            }
            for (int a13 = a12; a13 < n1[2]; a13++) {
              A1_pi0(2, _) = rep(pbinomial_pi0(a13, n1[2]), v);
              A1_pi1(2, _) = rep(pbinomial_pi1(a13, n1[2]), v);
              if (sum(A1_pi1) >= beta*pow(v, 2)) {
                break;
              }
              NumericVector a1 = NumericVector::create(a11, a12, a13);
              for (int n21 = max(1, Nmin - n11 - v + 1); n21 <= Nmax - n11 - v; n21++) {
                IntegerVector n2 = seq(n21, n21 + v - 1);
                for (int a21 = a11 + 1; a21 < n1[0] + n2[0]; a21++) {
                  NumericMatrix A2_pi0(v, v), A2_pi1(v, v);
                  for (int i = 0; i <= v - 1; i++) {
                    for (int s1 = a1[i] + 1; s1 <= min(a21, n1[i]); s1++) {
                      A2_pi0(i, 0) += dbinomial_pi0(s1, n1[i])*pbinomial_pi0(a21 - s1, n2[0]);
                      A2_pi1(i, 0) += dbinomial_pi1(s1, n1[i])*pbinomial_pi1(a21 - s1, n2[0]);
                    }
                  }
                  if (sum(A1_pi1) + sum(A2_pi1) >= beta*pow(v, 2)) {
                    break;
                  }
                  for (int a22 = max(a12 + 1, a21); a22 < n1[0] + n2[1]; a22++) {
                    for (int i = 0; i <= v - 1; i++) {
                      A2_pi0(i, 1) = 0;
                      A2_pi1(i, 1) = 0;
                      for (int s1 = a1[i] + 1; s1 <= min(a22, n1[i]); s1++) {
                        A2_pi0(i, 1) += dbinomial_pi0(s1, n1[i])*pbinomial_pi0(a22 - s1, n2[1]);
                        A2_pi1(i, 1) += dbinomial_pi1(s1, n1[i])*pbinomial_pi1(a22 - s1, n2[1]);
                      }
                    }
                    if (sum(A1_pi1) + sum(A2_pi1(_, 0)) + sum(A2_pi1(_, 1)) >= beta*pow(v, 2)) {
                      break;
                    }
                    for (int a23 = max(a13 + 1, a22); a23 < n1[0] + n2[2]; a23++) {
                      for (int i = 0; i <= v - 1; i++) {
                        A2_pi0(i, 2) = 0;
                        A2_pi1(i, 2) = 0;
                        for (int s1 = a1[i] + 1; s1 <= min(a23, n1[i]); s1++) {
                          A2_pi0(i, 2) += dbinomial_pi0(s1, n1[i])*pbinomial_pi0(a23 - s1, n2[2]);
                          A2_pi1(i, 2) += dbinomial_pi1(s1, n1[i])*pbinomial_pi1(a23 - s1, n2[2]);
                        }
                      }
                      if (sum(A1_pi1) + sum(A2_pi1) >= beta*pow(v, 2)) {
                        break;
                      }
                      if (mean(A1_pi0) + mean(A2_pi0) > 1 - alpha) {
                        double SESS_pi0 = 0;
                        double SESS_pi1 = 0;
                        for (int i = 1; i <= v; i++) {
                          for (int j = 1; j <= v; j++) {
                            SESS_pi0 += n1[i - 1] + (1 - A1_pi0(i - 1, j - 1))*n2[j - 1];
                            SESS_pi1 += n1[i - 1] + (1 - A1_pi1(i - 1, j - 1))*n2[j - 1];
                          }
                        }
                        feasible_designs(counter, _) = NumericVector::create(n1[0], n1[1], n1[2],
                                         n2[0], n2[1], n2[2],
                                                  a11, a12, a13,
                                                  a21, a22, a23,
                                                  1 - mean(A1_pi0) - mean(A2_pi0),
                                                  1 - mean(A1_pi1) - mean(A2_pi1), SESS_pi0/pow(v, 2), SESS_pi1/pow(v, 2));
                        counter++;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                          counter - 1 : 0)),
                                                          Range(0, (J == 1 ? 2*v + 1 : 4*v + 3)));
  return output;
}



