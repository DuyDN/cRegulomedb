all: build_db clean compress_db upload_db
        
build_db:
	Rscript R/build_script.R
	Rscript R/db_testset.R
clean:
	rm -rf tmp

compress_db:
	gzip cRegulome.db
	gzip test.db

upload_db:
	~/dropbox_uploader.sh upload cRegulome.db.gz cRegulome.db.gz
	~/dropbox_uploader.sh upload test.db.gz test.db.gz
