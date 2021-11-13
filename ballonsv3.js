		var CLASS_AIDE = "aide";
		var CLASS_AIDESMALL = "aideSmall";
		var g_pos_curseur;
		var g_bulle_aide;
		var g_textes_bulle = new Array();
		var g_positions = new Array();
		var g_bulle_flag = false;
		var g_id_cour = null;
		onload = init;


    // Get the y position of an element. Needs to loop through all the parents
    function getYPos(el) 
    {
        for (var lx=0, ly=0;
             el != null;
             lx += el.offsetLeft, ly += el.offsetTop, el = el.offsetParent);
        return ly;
    } 

    // Get the child of type IMG
    function getIMGChild(el)
    {	
    	for (var i=0; i<el.childNodes.length; ++i)
		{	
			var x=el.childNodes[i];			
			if (x.tagName == "IMG")									
    			return el.childNodes[i];
    	}
    	return "";
	}

    

		function init()
		{

			// creation du <div> contenant la bulle d'aide
			g_bulle_aide = creer_bulle_aide();
			// affectation de onmouseover et de onmouseout sur tous les <span class="aide">
			var tab_spans = document.getElementsByTagName("span");
			for (i in tab_spans)
			{
				if (tab_spans[i].className == CLASS_AIDE || tab_spans[i].className == CLASS_AIDESMALL)
				{
					tab_spans[i].id = CLASS_AIDE + i;
					tab_spans[i].onmouseover = afficher_aide;
					// astuce pour le cas o� le title contient un lien :
					// seuls les textes ne contenant pas de lien ferment la bulle sur un onmouseout
					// pour les autres, on fait un 2e onmouseover
					if (!tab_spans[i].title.sansMouseOut())
						tab_spans[i].onmouseout = masquer_aide;
					// on stocke les title dans un tableau, et on supprime chaque title
					// pour ne pas les afficher en m�me temps que la bulle d'aide
					g_textes_bulle[CLASS_AIDE + i] = tab_spans[i].title;
					imgh = getIMGChild(tab_spans[i]);
					// imgh = tab_spans[i];
					imgh = imgh.clientHeight;
					// g_textes_bulle[CLASS_AIDE + i] = imgh  +"px";
					g_positions[CLASS_AIDE + i] = getYPos(tab_spans[i]) - imgh  +"px";
					// alert(g_positions[CLASS_AIDE + i]);
					tab_spans[i].title = "";
				}
			}
		}
		function creer_bulle_aide()
		{
			return document.getElementById("bulle_aide");

		}
		function afficher_aide()
		{
			// affiche l'aide dans la bulle en fonction de l'id du span
			var texte = g_textes_bulle[this.id];
			var span_position = g_positions[this.id];

			// si le flag est true, on masquera la bulle (2e mouseover) si le title du span contient un lien
			if (g_bulle_flag)
			{
				if (this.id == g_id_cour)
				{
					var afficher = false;
					g_bulle_flag = false;
					g_id_cour = null;
				}
				else
				{
					var afficher = true;
					if (texte.sansMouseOut()) g_id_cour = this.id;
					else
					{
						g_bulle_flag = false;
						g_id_cour = null;
					}
				}
			}
			// sinon, on affichera la bulle
			// si le title du span contient un lien : on met le flag � true (1e mouseover)
			else
			{
				var afficher = true;
				if (texte.sansMouseOut())
				{
					g_bulle_flag = true;
					g_id_cour = this.id;
				}
			}
			// si affichage de la bulle
			if (afficher)
			{
				g_bulle_aide.innerHTML = texte;
				g_bulle_aide.style.visibility = "visible";
				g_bulle_aide.style.width = "50%";
				g_bulle_aide.style.left = "300px";
				g_bulle_aide.style.top = span_position;
				// g_bulle_aide.style.top = mouseY;
				this.style.color = "#990000";
			}
			// sinon masquage de la bulle
			else
				g_bulle_aide.style.visibility = "hidden";

		}
		function masquer_aide()
		{
			// masque la bulle d'aide (d�clench� par un onmouseout sur le <span>)
			g_bulle_aide.style.visibility = "hidden";
			this.style.color = "#000000";

		}

		String.prototype.sansMouseOut = function()
		{
			// d�termine si le texte de la bulle contient certaines balises HTML (actuellement seulement des liens)
			// qui emp�chent la fermeture de la bulle sur un onmousout du span ayant ouvert la bulle
			if (this.search(/<a href=[^>]+>.*<\/a>/) == -1)
				return false;
			else
				return true;
		}

		// affectation de la fonction "emplacementSouris" � l'�v�nement onmousemove
		document.onmousemove = emplacementSouris;

		// Cette fonction appelle emplacementSouris � chaque mouvement de la souris
		function emplacementSouris(e)
		{
			// avec IE, on utilise l'objet "event"
			if (document.all) g_pos_curseur = {x: event.offsetX, y: event.offsetY};
			// avec Netscape, Mozilla, on utilise l'�v�nement "e" en argument
			else g_pos_curseur = {x: e.pageX, y: e.pageY};
		}
