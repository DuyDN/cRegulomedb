all: build_db

build_db:
	R CMD BATCH R/build_mir.R
	R CMD BATCH R/build_tf.R
	R CMD BATCH R/build_testset.R

compress_db:
	gzip cRegulome.db
	gzip test.db
