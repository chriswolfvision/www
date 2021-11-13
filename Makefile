all:
	@echo "usage: make publish"

publish:
	rsync -vrout --delete --progress --modify-window=2 -l --exclude=.mypasswds --exclude=.htaccess --exclude=.svn-access-file --exclude=.ssh . cwolf@connect.liris.cnrs.fr:.
	#rsync -vrout --delete --progress --modify-window=2 -l --exclude=.mypasswds --exclude=.htaccess --exclude=.svn-access-file --exclude=.ssh . cwolf@liris.cnrs.fr:.

check:
	ls papers/* | ./scripts/findunref.sh .
	ls miscphotos/* | ./scripts/findunref.sh .
	ls graphics/* | ./scripts/findunref.sh .
