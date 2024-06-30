all:
	@echo "usage: make publish"
	@echo "usage: make clearcredentials"

publish:
	git commit -a -m "stuff"
	git push

check:
	ls papers/* | ./scripts/findunref.sh .
	ls miscphotos/* | ./scripts/findunref.sh .
	ls graphics/* | ./scripts/findunref.sh .

clearcredentials:
	git config credential.helper store
