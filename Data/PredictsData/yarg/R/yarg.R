.onLoad <- function(libname, pkgname) {
    # DBExtractDate and DBExtractDateTime are used by ReadDBExtract()
    # http://stackoverflow.com/questions/13022299/specify-date-format-for-colclasses-argument-in-read-table-read-csv
    setClass('DBExtractDate')
    setAs("character", "DBExtractDate", function(from) as.Date(from, format='%Y-%m-%d'))

    setClass('DBExtractDateTime')
    setAs("character", "DBExtractDateTime", function(from) as.POSIXct(from, format='%Y-%m-%d %H:%M:%S'))
}
