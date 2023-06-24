all:
	@echo "usage: make publish"

publish:
	git commit -a -m "stuff"
	git push

check:
	ls papers/* | ./scripts/findunref.sh .
	ls miscphotos/* | ./scripts/findunref.sh .
	ls graphics/* | ./scripts/findunref.sh .
