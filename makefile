all: build_db
        
build_db:
	Rscript R/build_script.R

clean:
	rm -rf tmp

