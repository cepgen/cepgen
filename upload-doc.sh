./bin/cepgenDocGenerator -o ../doc/.raw_modules.html -e
ninja Sphinx
rsync -arvz -e 'ssh -p 222' --delete doc/output_html/ lforthomme@login.hepforge.org:~/cepgen/public_html/
