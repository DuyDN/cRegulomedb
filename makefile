all: build_db clean compress_db upload_db
        
build_db:
	Rscript R/build_script.R

clean:
	rm -rf tmp

compress_db:
	gzip cRegulome.db.gz

upload_db:
	~/dropbox_uploader.sh chunck_upload cRegulome.db.gz cRegulome.db.gz
