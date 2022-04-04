./cepgenDocGenerator -o doc/output_html/.raw_modules.html -b 0
rsync -arvz -e 'ssh -p 222' --delete doc/output_html/ lforthomme@login.hepforge.org:~/cepgen/public_html/
