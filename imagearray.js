var ARR = new Array();
ARR[0] = new Image(); ARR[0].src = "miscphotos/chris_usa_yosemite_2014.jpg";
ARR[1] = new Image(); ARR[1].src = "miscphotos/chris-zr7-silhouette.jpg";
ARR[2] = new Image(); ARR[2].src = "miscphotos/chris_chateau_2012_corr.jpg";
ARR[3] = new Image(); ARR[3].src = "miscphotos/chris_toulon2015.jpg";
ARR[4] = new Image(); ARR[4].src = "miscphotos/chris-noah-rue-merciere.jpg";
ARR[5] = new Image(); ARR[5].src = "miscphotos/chris-icpr2012-www.jpg";
ARR[6] = new Image(); ARR[6].src = "miscphotos/chris_bresse_120723.jpg";
ARR[7] = new Image(); ARR[7].src = "miscphotos/chris_noah_salzburg.jpg";
ARR[8] = new Image(); ARR[8].src = "miscphotos/chris3.jpg";
ARR[9] = new Image(); ARR[9].src = "miscphotos/chris5.jpg";
ARR[10] = new Image(); ARR[10].src = "miscphotos/chris7.jpg";
ARR[11] = new Image(); ARR[11].src = "miscphotos/chris11.jpg";
ARR[12] = new Image(); ARR[12].src = "miscphotos/chris12.jpg";
ARR[13] = new Image(); ARR[13].src = "miscphotos/chris15.jpg";
curIndex=0;

function nextImageInArray() 
{
	curIndex=curIndex+1;
	if (curIndex>13)
		curIndex=0;
	node = document.getElementById("imagearraytarget");
	node.src=ARR[curIndex].src;
}

		
