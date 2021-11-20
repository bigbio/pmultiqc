$(document).ready(function () {
    var quant_table = document.getElementById("quantification_of_peptides");
    var quantTable = quant_table.getElementsByTagName("tbody")[0];
    var quantTotalPage = document.getElementById("quantTotalPage");
    var quantPageNum = document.getElementById("quantPageNum");

    quantPreDom = document.getElementById("quantPre");
    quantNextDom = document.getElementById("quantNext");
    quantFirstDom = document.getElementById("quantFirst");
    quantLastDom = document.getElementById("quantLast");
    quant_numrows = document.getElementById("quantification_of_peptides_numrows_text");
    quant_sub = quant_numrows.getElementsByTagName("sub")[0];
    quant_trs = quantTable.getElementsByTagName("tr");

    numberRowsInQuantTable = 0;
    quantPageSize = 50;
    quantPage = 1;
    quantPageCount().then(res => { quantTotalPage.innerHTML = parseInt(res / 50 + 1);
        quant_sub.innerHTML = res;
        numberRowsInQuantTable = res;
        quantLastRows = quantPageSize * (parseInt(res / 50 + 1) - 1);
    });
    quantNextLink();
    quantLastLink();
    quantPageNum.innerHTML = '1';
})

//NextPage
function quantNext(){
    currentRow = quantPageSize * quantPage;
    maxRow = currentRow + quantPageSize;
    if ( maxRow > numberRowsInQuantTable ) maxRow = numberRowsInQuantTable;
    updateQuantData(currentRow).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = quant_trs[i].getElementsByTagName("td");
			quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
	});
    quantPage++;

    if ( maxRow == numberRowsInQuantTable ) { quantNextText(); quantLastText(); }
    showPage(quantPageNum, quantPage);
    quantPreLink();
    quantFirstLink();
}

//PreviousPage
function quantPre(){
    quantPage--;
    currentRow = quantPageSize * quantPage;
    maxRow = currentRow - quantPageSize;
    if ( currentRow > numberRowsInQuantTable ) currentRow = numberRowsInQuantTable;
	updateQuantData(currentRow - quantPageSize).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = quant_trs[i].getElementsByTagName("td");
			quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
	});


    if ( maxRow === 0 ){ quantPreText(); quantFirstText(); }
    showPage(quantPageNum, quantPage);
    quantNextLink();
    quantLastLink();
}

//FirstPage
function quantFirst(){
    quantPage = 1;
    updateQuantData(0).then(res =>{
    	for(i=0; i<res.length; i++){
    		console.log(res[i]);
    		tds = quant_trs[i].getElementsByTagName("td");
    		quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
    		for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
    		}
    	}
		console.log(i);
		for (k = 1; k < 50; k++){
			quant_trs[k].style.display = '';
		}

    });
    showPage(quantPageNum, quantPage);

    quantPreText();
    quantNextLink();
    quantLastLink();
}

//LastPage
function quantLast(){
    quantPage = parseInt(quantLastRows / quantPageSize + 1);
    updateQuantData(quantLastRows).then(res =>{
    	for(i=0; i<res.length; i++){
    		console.log(res[i]);
    		tds = quant_trs[i].getElementsByTagName("td");
    		quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
    		for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
    		}
    	}
		if ( i <= quant_trs.length - 1){
			for (k = i; k < quant_trs.length; k++){
				        quant_trs[k].style.display = 'none';
			}
		}
    });
    showPage(quantPageNum, quantPage);

    quantPreLink();
    quantNextText();
    quantFirstLink();
}

function showPage(pageNum, page){
    pageNum.innerHTML = page;

}

//TotalPage
async function quantPageCount(){
    let t;
    await axios.get("proteomicslfq.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
        .then(function (response) {
        let db = new window.SQL.Database(new Uint8Array(response.data));
        // execute query
        let s = new Date().getTime();
        let r = db.exec("select count(*) from quant");
        let e = new Date().getTime();
        console.info("Time consuming to query data：" + (e - s) + "ms");
        // parse data
        console.info(r[0]['values'][0][0]);
        t = r[0]['values'][0][0];
        })
        .catch(function (error) {
            console.info(error);
    });
    return t;
}

//TotalPage
async function updateQuantData(currentRow){
	let d;
    await axios.get("proteomicslfq.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*',
			'Access-Control-Allow-Headers': '*'}})
    		.then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));
			// 执行查询
			let s = new Date().getTime();
			let r = db.exec("select * from quant " + "limit "+ String(currentRow) + ",50");
			let e = new Date().getTime();
			console.info("Time consuming to query data：" + (e - s) + "ms");
			// parse data
			console.info(r);
			d = r[0]['values'];
    })
    .catch(function (error) {
                console.info(error);
    });
	return d;
}

//ShowLink
function quantPreLink(){ quantPreDom.innerHTML = "<a href='javascript:quantPre();'>Previous Page</a>";}
function quantPreText(){ quantPreDom.innerHTML = "Previous Page";}

function quantNextLink(){ quantNextDom.innerHTML = "<a href='javascript:quantNext();'>Next Page</a>";}
function quantNextText(){ quantNextDom.innerHTML = "Next Page";}

function quantFirstLink(){ quantFirstDom.innerHTML = "<a href='javascript:quantFirst();'>First Page</a>";}
function quantFirstText(){ quantFirstDom.innerHTML = "First Page";}

function quantLastLink(){ quantLastDom.innerHTML = "<a href='javascript:quantLast();'>Last Page</a>";}
function quantLastText(){ quantLastDom.innerHTML = "Last Page";}

function quant_page_jump(){
	if(event.keyCode === 13){
		quantPage = document.getElementById("pep_page").value;
		if (quantPage > parseInt(numberRowsInQuantTable / 50 + 1) || quantPage === ""){
			alert("not valid page!");
		}

		else if(quantPage === parseInt(numberRowsInQuantTable / 50 + 1)){
			quantLast();
		} else{
			currentRow = quantPageSize * quantPage;
			maxRow = currentRow - quantPageSize;
			updateQuantData(currentRow - quantPageSize).then(res =>{
  				for(i=0; i<res.length; i++){
  					console.log(res[i]);
  					tds = quant_trs[i].getElementsByTagName("td");
  					quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
  					for(j=0; j<tds.length; j++){
  						if(res[i][j + 1] == null){
					        tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
  						}
				        else{
					        tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				        }
  					}
  				}
  				console.log(i);
  				for (k = 1; k < 50; k++){
  					quant_trs[k].style.display = '';
  				}

  			});

			if ( maxRow === 0 ){ quantPreText(); quantFirstText(); }
			showPage(quantPageNum, quantPage);
			quantPreLink();
			quantNextLink();
			quantFirstLink();
			quantLastLink();
		}
	}
}

async function searchData(filter, col, table){
	let d;
    await axios.get("proteomicslfq.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
    		.then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));

			// execute query
			let s = new Date().getTime();
			let r = db.exec("select * from "+ table  +" where "+ col + " like '%" + String(filter) + "%'");
			let e = new Date().getTime();
			console.info("Time consuming to query data：" + (e - s) + "ms");
			// parse data
            if (r.length == 0) {d = r;}
			else{d = r[0]['values'];}
			console.log(d);
    })
    .catch(function (error) {
        console.info(error);
    });
	return d;
}




function searchQuantFunction() {
    if (event.keyCode === 13) {
        var myInput=document.getElementById("quant_search");
        var filter=myInput.value.toUpperCase();
        var search_col=document.getElementById("quant_search_col");
        var index = search_col.selectedIndex;
        var value = search_col.options[index].text;

        searchData(filter, value, 'quant').then(res =>{
        for(i=0; i<res.length; i++){
            if(i>=50){
                break;
            }
            quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            tds = quant_trs[i].getElementsByTagName("td");
            for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
                    tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                }
                else{
                    tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                }
            }

        }
        for (k = 0; k < 50; k++){
            if (k>=i) {
                quant_trs[k].style.display = 'none';
            } else{
                quant_trs[k].style.display = '';
            }
            }
        });
    }
}




$(document).ready(function () {
    var psm_table = document.getElementById("peptide_spectrum_match");
    psmTable = psm_table.getElementsByTagName("tbody")[0];

    psmTotalPage = document.getElementById("psmTotalPage");
    psmPageNum = document.getElementById("psmPageNum");

    psmPreDom = document.getElementById("psmPre");
    psmNextDom = document.getElementById("psmNext");
    psmFirstDom = document.getElementById("psmFirst");
    psmLastDom = document.getElementById("psmLast");
    psm_numrows = document.getElementById("peptide_spectrum_match_numrows_text");
    psm_sub = psm_numrows.getElementsByTagName("sub")[0];
    psm_trs = psmTable.getElementsByTagName("tr");

    numberRowsInPsmTable = 0;
    psmPageSize = 50;
    psmPage = 1;

    psmPageCount().then(res => { psmTotalPage.innerHTML = parseInt(res / 50 + 1);
        psm_sub.innerHTML = res;
        numberRowsInPsmTable = res;
        psmLastRows = psmPageSize * (parseInt(res / 50 + 1) - 1);
    });
    psmNextLink();
    psmLastLink();
    psmPageNum.innerHTML = '1';
})

    //NextPage
function psmNext(){
    currentRow = psmPageSize * psmPage;
    maxRow = currentRow + psmPageSize;
    if ( maxRow > numberRowsInPsmTable ) maxRow = numberRowsInPsmTable;
        updatePsmData(currentRow).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = psm_trs[i].getElementsByTagName("td");
			psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
	});
    psmPage++;

    if ( maxRow == numberRowsInPsmTable ) { psmNextText(); psmLastText(); }
    showPage(psmPageNum, psmPage);
    psmPreLink();
    psmFirstLink();
}

//PreviousPage
function psmPre(){
    psmPage--;
    currentRow = psmPageSize * psmPage;
    maxRow = currentRow - psmPageSize;
    if ( currentRow > numberRowsInPsmTable ) currentRow = numberRowsInPsmTable;
	updatePsmData(currentRow - psmPageSize).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = psm_trs[i].getElementsByTagName("td");
			psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
	});


    if ( maxRow === 0 ){ psmPreText(); psmFirstText(); }
    showPage(psmPageNum, psmPage);
    psmNextLink();
    psmLastLink();
}

//FirstPage
function psmFirst(){
    psmPage = 1;
    updatePsmData(0).then(res =>{
    	for(i=0; i<res.length; i++){
    		console.log(res[i]);
    		tds = psm_trs[i].getElementsByTagName("td");
    		psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
    		for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
    		}
    	}
		console.log(i);
		for (k = 1; k < 50; k++){
			quant_trs[k].style.display = '';
		}

    });
    showPage(psmPageNum, psmPage);

    psmPreText();
    psmNextLink();
    psmLastLink();
}

//LastPage
function psmLast(){
    psmPage = parseInt(psmLastRows / psmPageSize + 1);
    updatePsmData(psmLastRows).then(res =>{
    	for(i=0; i<res.length; i++){
    		console.log(res[i]);
    		tds = psm_trs[i].getElementsByTagName("td");
    		psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
    		for(j=0; j<tds.length; j++){
    			if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
    		}
    	}
		if ( i <= psm_trs.length - 1){
			for (k = i; k < psm_trs.length; k++){
				        psm_trs[k].style.display = 'none';
			}
		}
    });
    showPage(psmPageNum, psmPage);

    psmPreLink();
    psmNextText();
    psmFirstLink();
}

//TotalPage
async function psmPageCount(){
    let t;
    await axios.get("proteomicslfq.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
        .then(function (response) {
        let db = new window.SQL.Database(new Uint8Array(response.data));
        let r = db.exec("select count(*) from psm");
        t = r[0]['values'][0][0];
        })
        .catch(function (error) {
            console.info(error);
    });
    return t;
}

//ShowLink
function psmPreLink(){ psmPreDom.innerHTML = "<a href='javascript:psmPre();'>Previous Page</a>";}
function psmPreText(){ psmPreDom.innerHTML = "Previous Page";}

function psmNextLink(){ psmNextDom.innerHTML = "<a href='javascript:psmNext();'>Next Page</a>";}
function psmNextText(){ psmNextDom.innerHTML = "Next Page";}

function psmFirstLink(){ psmFirstDom.innerHTML = "<a href='javascript:psmFirst();'>First Page</a>";}
function psmFirstText(){ psmFirstDom.innerHTML = "First Page";}

function psmLastLink(){ psmLastDom.innerHTML = "<a href='javascript:psmLast();'>Last Page</a>";}
function psmLastText(){ psmLastDom.innerHTML = "Last Page";}

async function updatePsmData(currentRow){
	let d;
    await axios.get("proteomicslfq.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
    		.then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));
			let r = db.exec("select * from psm " + "limit "+ String(currentRow) + ",50");
			d = r[0]['values'];
    })
    .catch(function (error) {
        console.info(error);
    });
	return d;
}

function psm_page_jump(){
	if(event.keyCode === 13){
		psmPage = document.getElementById("psm_page").value;
		if (psmPage > parseInt(numberRowsInPsmTable / 50 + 1) || psmPage === ""){
			alert("not valid page!");
		}

		else if(psmPage === parseInt(numberRowsInPsmTable / 50 + 1)){
			psmLast();
		} else{
			currentRow = psmPageSize * psmPage;
			maxRow = currentRow - psmPageSize;
			updatePsmData(currentRow - psmPageSize).then(res =>{
  				for(i=0; i<res.length; i++){
  					tds = psm_trs[i].getElementsByTagName("td");
  					psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
  					for(j=0; j<tds.length; j++){
  						if(res[i][j + 1] == null){
					        tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
  						}
				        else{
					        tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				        }
  					}
  				}
  				for (k = 1; k < 50; k++){
  					psm_trs[k].style.display = '';
  				}

  			});

			if ( maxRow === 0 ){ psmPreText(); psmFirstText(); }
			showPage(psmPageNum, psmPage);
			psmPreLink();
			psmNextLink();
			psmFirstLink();
			psmLastLink();
		}
	}
}

function searchPsmFunction() {
    if (event.keyCode === 13) {
        var myInput=document.getElementById("psm_search");
        var filter=myInput.value.toUpperCase();
        var search_col=document.getElementById("psm_search_col");
        var index = search_col.selectedIndex;
        var value = search_col.options[index].text;

        searchData(filter, value, 'psm').then(res =>{
        for(i=0; i<res.length; i++){
            if(i>=50){
                break;
            }
            psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            tds = psm_trs[i].getElementsByTagName("td");
            for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
                    tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                }
                else{
                    tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                }
            }

        }
        for (k = 0; k < 50; k++){
            if (k>=i) {
                psm_trs[k].style.display = 'none';
            } else{
                psm_trs[k].style.display = '';
            }
            }
        });
    }
}