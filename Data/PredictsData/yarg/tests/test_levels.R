TestEnsureFactorLevels <- function() {
  AssertEqual(LETTERS[1:4], levels(EnsureFactorLevels(LETTERS[1], LETTERS[1:4])))
  AssertEqual(LETTERS[1:4], levels(EnsureFactorLevels(LETTERS[2:4], LETTERS[1:4])))
  AssertEqual(LETTERS[1:4], levels(EnsureFactorLevels(LETTERS[1:4], LETTERS[1:4])))
  AssertRaises(EnsureFactorLevels(LETTERS[1:5], LETTERS[1:4]))
}
