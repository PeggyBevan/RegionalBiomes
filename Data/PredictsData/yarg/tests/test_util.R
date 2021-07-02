TestParseRangeOrValue <- function() {
    AssertEqual(c(NA,NA,1),as.vector(ParseRangeOrValue('1')))
    AssertEqual(as.numeric(c(NA,NA,NA)),as.vector(ParseRangeOrValue('')))
    AssertEqual(c(1,2,NA),as.vector(ParseRangeOrValue('1-2')))
    AssertEqual(c(10,20,NA),as.vector(ParseRangeOrValue('  10 - 20  ')))

    AssertRaises(ParseRangeOrValue('?', print.errors=FALSE))
    AssertRaises(ParseRangeOrValue('1.1.1', print.errors=FALSE))
    AssertRaises(ParseRangeOrValue('-1', print.errors=FALSE))
    AssertRaises(ParseRangeOrValue('1-a', print.errors=FALSE))
    AssertRaises(ParseRangeOrValue('a', print.errors=FALSE))
}

TestMeanOfRangeOrValue <- function() {
    AssertEqual(1, MeanOfRangeOrValue('1'))
    AssertEqual(as.numeric(NA), MeanOfRangeOrValue(''))
    AssertEqual(1.5, MeanOfRangeOrValue('1-2'))
    AssertEqual(15, MeanOfRangeOrValue('  10 - 20  '))
}

TestDiffOfRangeOrValue <- function() {
    AssertEqual(1, DiffOfRangeOrValue('1'))
    AssertEqual(as.numeric(NA), DiffOfRangeOrValue(''))
    AssertEqual(1, DiffOfRangeOrValue('1-2'))
    AssertEqual(10, DiffOfRangeOrValue('  10 - 20  '))
}
