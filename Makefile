all:
	@echo "usage: make publish"

publish:
	echo "do a commit and push!"

check:
	ls papers/* | ./scripts/findunref.sh .
	ls miscphotos/* | ./scripts/findunref.sh .
	ls graphics/* | ./scripts/findunref.sh .
