// --------------------------------------------------------------------------
// The scripts which display the abstract or the BibTeX entry of a 
// publication on a mouse click.
// --------------------------------------------------------------------------

function displayElementOf(id,element) 
{
	// Get the main node
	mnode=document.getElementById(id);
	
	// Travers the children and search for the abstract,
	// i.e. the node with class "t-abstract"
	for (var i=0; i<mnode.childNodes.length; ++i)
	{
		var x=mnode.childNodes[i];			
		if (x.className == element)									
			if (x.style.display != 'block')
				x.style.display='block';
			else
				x.style.display='none';
	}		
}

function displayAbstractOf(id)
{
	displayElementOf(id,'t-abstract');
}

function displayBibtexOf(id)
{
	displayElementOf(id,'t-bibtex');
}

// --------------------------------------------------------------------------
// marks an object with two red bars
// --------------------------------------------------------------------------

function markObj(id) 
{
	identity=document.getElementById(id);
	identity.style.borderLeft='15px solid #880000';
	identity.style.borderRight='15px solid #880000';
	identity.style.paddingLeft='2em';
	identity.style.paddingRight='2em';
}

// --------------------------------------------------------------------------
// The script which is executed at each page loading 
// --------------------------------------------------------------------------

function startPage() {
	xpos=location.href.indexOf("#");
	if (xpos!=-1) 
	{
		id=location.href.substr(xpos+1);
		displayAbstractOf(id);		
		markObj(id);			
		
		/* Move to the object position */
		abs=document.getElementById(id);
		window.scrollTo(abs.offsetLeft,abs.offsetTop);
	}
}
