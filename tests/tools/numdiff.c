/*
 * numdiff -- compare two text files numerically.
 *
 *   numdiff <reference> <output> [rel_tol] [abs_tol]
 *
 * Files are split into tokens on whitespace and commas. Tokens that parse
 * as numbers must satisfy |a - b| <= abs_tol + rel_tol * max(|a|, |b|);
 * other tokens must match exactly. Exits 0 if the files agree, 1 if not,
 * 2 on usage/IO errors. Reports the first few differences.
 */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int next_token(FILE *f, char *buf, size_t n, long *line) {
  int c;
  size_t k = 0;
  while ((c = fgetc(f)) != EOF && (isspace(c) || c == ','))
    if (c == '\n') ++*line;
  if (c == EOF) return 0;
  do {
    if (k + 1 < n) buf[k++] = (char)c;
  } while ((c = fgetc(f)) != EOF && !isspace(c) && c != ',');
  if (c == '\n') ungetc(c, f);
  buf[k] = '\0';
  return 1;
}

static int as_number(const char *s, double *v) {
  char *end;
  *v = strtod(s, &end);
  return end != s && *end == '\0';
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s <reference> <output> [rel_tol] [abs_tol]\n", argv[0]);
    return 2;
  }
  double rel = argc > 3 ? atof(argv[3]) : 1e-8;
  double abs_tol = argc > 4 ? atof(argv[4]) : 1e-12;
  FILE *a = fopen(argv[1], "r"), *b = fopen(argv[2], "r");
  if (!a || !b) {
    fprintf(stderr, "numdiff: cannot open %s\n", !a ? argv[1] : argv[2]);
    return 2;
  }
  char ta[256], tb[256];
  long la = 1, lb = 1, ntok = 0, ndiff = 0;
  double worst = 0;
  for (;;) {
    int ga = next_token(a, ta, sizeof ta, &la), gb = next_token(b, tb, sizeof tb, &lb);
    if (!ga && !gb) break;
    if (ga != gb) {
      fprintf(stderr, "numdiff: %s has more tokens (line %ld)\n", ga ? argv[1] : argv[2], ga ? la : lb);
      return 1;
    }
    ++ntok;
    double x, y;
    int ok;
    if (as_number(ta, &x) && as_number(tb, &y)) {
      double d = fabs(x - y), tol = abs_tol + rel * fmax(fabs(x), fabs(y));
      ok = d <= tol || (isnan(x) && isnan(y));
      if (fmax(fabs(x), fabs(y)) > 0) worst = fmax(worst, d / fmax(fabs(x), fabs(y)));
    } else {
      ok = strcmp(ta, tb) == 0;
    }
    if (!ok && ++ndiff <= 5)
      fprintf(stderr, "numdiff: line %ld: reference '%s' vs output '%s'\n", la, ta, tb);
  }
  printf("numdiff: %ld tokens compared, %ld differ, max relative difference %.3g\n", ntok, ndiff, worst);
  return ndiff ? 1 : 0;
}
