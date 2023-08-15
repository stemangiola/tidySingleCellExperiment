df <- pbmc_small
df$number <- sample(seq(ncol(df)))
df$factor <- sample(gl(2, 1, ncol(df), c("g1", "g2")))

test_that("un/nest()", {
    fd <- nest(df, data=-factor)
    expect_equal(
        nrow(fd), 
        nlevels(df$factor))
    expect_equal(
        vapply(fd$data, ncol, integer(1)), 
        tabulate(df$factor))
    fd <- unnest(fd, data)
    # [HLC: this is necessary because unnest() 
    # currently duplicates the 'int_metadata']
    int_metadata(fd) <- int_metadata(df)
    expect_equal(fd, df[, colnames(fd)])
})

test_that("unite()/separate()", {
    expect_error(unite(df, "x", c(number, x)))
    expect_error(separate(df, x, c("a", "b")))
    
    fd <- unite(df, "string", c(number, factor), sep=":")
    expect_null(fd$number)
    expect_null(fd$factor)
    expect_identical(fd$string, paste(df$number, df$factor, sep=":"))
    
    fd <- separate(fd, string, c("a", "b"), sep=":")
    expect_null(fd$string)
    expect_identical(fd$a, paste(df$number))
    expect_identical(fd$b, paste(df$factor))
    
    # special columns are blocked
    expect_error(unite(df, ".cell", c(number, factor), sep=":"))
    fd <- df; colnames(fd) <- paste(colnames(df), "x", sep="-")
    expect_error(separate(fd, .cell, c("a", "b"), sep="-"))
})

test_that("extract()", {
    expect_error(extract(df, a, "b"))
    expect_error(extract(df, factor, "x", ""))
    
    fd <- mutate(df, factor=paste(factor))
    expect_identical(extract(df, factor, "factor"), fd)
    expect_identical(extract(df, factor, "factor", "(.*)"), fd)
    
    fd <- extract(df, factor, "g", "g([0-9])", convert=FALSE)
    expect_identical(fd$g, gsub("^g", "", df$factor))
    expect_null(fd$factor)
    
    fd <- extract(df, factor, "g", "g([0-9])", convert=TRUE)
    expect_identical(fd$g, as.integer(gsub("^g", "", df$factor)))
    expect_null(fd$factor)
})

test_that("pivot_longer()", {
    abc <- c("a", "b", "c")
    df$string <- sample(abc, ncol(df), TRUE)
    fd <- pivot_longer(df, c(factor, string))
    expect_s3_class(fd, "tbl_df")
    expect_true(!any(c("factor", "string") %in% names(fd)))
    expect_setequal(unique(fd$name), c("factor", "string"))
    expect_setequal(unique(fd$value), c(levels(df$factor), abc))
})
