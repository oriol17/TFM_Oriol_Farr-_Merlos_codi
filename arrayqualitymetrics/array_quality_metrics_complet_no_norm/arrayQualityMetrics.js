// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, true, true ];
var arrayMetadata    = [ [ "1", "GSM1134016_EA04058_25073_H133+_F15.CEL", "1" ], [ "2", "GSM1134017_EA04058_25074_H133+_F16.CEL", "2" ], [ "3", "GSM1134018_EA04058_25075_H133+_F20.CEL", "3" ], [ "4", "GSM1134019_EA04058_25076_H133+_F27.CEL", "4" ], [ "5", "GSM1134020_EA04058_25077_H133+_F34.CEL", "5" ], [ "6", "GSM1134021_EA04058_25078_H133+_F38.CEL", "6" ], [ "7", "GSM1134022_EA04058_25079_H133+_F42.CEL", "7" ], [ "8", "GSM1134023_EA04058_25080_H133+_F43.CEL", "8" ], [ "9", "GSM1134024_EA04058_25081_H133+_F45.CEL", "9" ], [ "10", "GSM1134025_EA04058_25082_H133+_F52.CEL", "10" ], [ "11", "GSM1134026_EA04058_28124_H133+_TS208.CEL", "11" ], [ "12", "GSM1134027_EA04058_28335_H133+_TS250.CEL", "12" ], [ "13", "GSM1134028_EA04058_30434_H133+_TS258.CEL", "13" ], [ "14", "GSM1134029_EA04058_30436_H133+_TS245.CEL", "14" ], [ "15", "GSM1134030_EA04058_30437_H133+_TS238.CEL", "15" ], [ "16", "GSM1134031_EA04058_52889_H133+_TS255.CEL", "16" ], [ "17", "GSM1134032_EA04058_52892_H133+_TS279.CEL", "17" ], [ "18", "GSM1134033_EA04058_52893_H133+_TS236.CEL", "18" ], [ "19", "GSM1134034_EA04058_52894_H133+_TS285.CEL", "19" ], [ "20", "GSM1134035_EA04058_52895_H133+_TS240.CEL", "20" ], [ "21", "GSM1134036_EA04058_52896_H133+_TS299.CEL", "21" ], [ "22", "GSM1134037_EA04058_52898_H133+_TS108.CEL", "22" ], [ "23", "GSM1134038_EA04058_52899_H133+_TS221.CEL", "23" ], [ "24", "GSM1134039_EA04058_53536_H133+_TS13.CEL", "24" ], [ "25", "GSM1134040_EA04058_53538_H133+_TS290.CEL", "25" ], [ "26", "GSM1134041_EA04058_53540_H133+_TS181.CEL", "26" ], [ "27", "GSM1134042_EA04058_30426_H133+_TS168.CEL", "27" ], [ "28", "GSM1134043_EA04058_30427_H133+_TS79.CEL", "28" ], [ "29", "GSM1134044_EA04058_52884_H133+_TS11.CEL", "29" ], [ "30", "GSM1134045_EA04058_52885_H133+_TS120.CEL", "30" ], [ "31", "GSM1134046_EA04058_52886_H133+_TS189.CEL", "31" ], [ "32", "GSM1134047_EA04058_52887_H133+_TS234.CEL", "32" ], [ "33", "GSM1134048_EA04058_52888_H133+_TS241.CEL", "33" ], [ "34", "GSM1134049_EA04058_52890_H133+_TS266.CEL", "34" ], [ "35", "GSM1134050_EA04058_52891_H133+_TS315.CEL", "35" ], [ "36", "GSM1134051_EA04058_52897_H133+_TS308.CEL", "36" ], [ "37", "GSM1411021_Turner_1_HG-U133_Plus_2_.CEL", "37" ], [ "38", "GSM1411022_Turner_2_HG-U133_Plus_2_.CEL", "38" ], [ "39", "GSM1411023_Turner_10_B_HG-U133_Plus_2_.CEL", "39" ], [ "40", "GSM1411024_Turner_12_B_HG-U133_Plus_2_.CEL", "40" ], [ "41", "GSM1411025_Turner_15_HG-U133_Plus_2_.CEL", "41" ], [ "42", "GSM1411026_Turner_11_control_HG-U133_Plus_2_.CEL", "42" ], [ "43", "GSM1411027_Turner_7_control_HG-U133_Plus_2_.CEL", "43" ], [ "44", "GSM1411028_Turner_13_B_control_HG-U133_Plus_2_.CEL", "44" ], [ "45", "GSM1411029_Turner_14_B_control_HG-U133_Plus_2_.CEL", "45" ], [ "46", "GSM1411030_Turner_9_control_HG-U133_Plus_2_.CEL", "46" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
