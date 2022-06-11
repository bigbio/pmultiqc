function draw_sparkline(){
    Highcharts.SparkLine = function (a, b, c) {
        const hasRenderToArg = typeof a === 'string' || a.nodeName;
        let options = arguments[hasRenderToArg ? 1 : 0];
        const defaultOptions = {
            chart: {
            renderTo: (options.chart && options.chart.renderTo) || (hasRenderToArg && a),
            backgroundColor: null,
            borderWidth: 0,
            type: 'area',
            margin: [2, 0, 2, 0],
            width: 120,
            height: 20,
            style: {
                overflow: 'visible'
            },
            // small optimalization, saves 1-2 ms each sparkline
            skipClone: true
            },
            title: {
            text: ''
            },
            credits: {
            enabled: false
            },
            xAxis: {
            labels: {
                enabled: false
            },
            title: {
                text: null
            },
            startOnTick: false,
            endOnTick: false,
            tickPositions: []
            },
            yAxis: {
            endOnTick: false,
            startOnTick: false,
            labels: {
                enabled: false
            },
            title: {
                text: null
            },
            tickPositions: [0]
            },
            legend: {
            enabled: false
            },
            tooltip: {
            hideDelay: 0,
            outside: true,
            shared: true
            },
            exporting: {
                enabled: false
            },
            plotOptions: {
            series: {
                animation: false,
                lineWidth: 1,
                shadow: false,
                states: {
                hover: {
                    lineWidth: 1
                }
                },
                marker: {
                radius: 1,
                states: {
                    hover: {
                    radius: 2
                    }
                }
                },
                fillOpacity: 0.25
            },
            column: {
                negativeColor: 'rgb(163,237,186)',
                borderColor: 'silver'
            }
            }
        };
        
        options = Highcharts.merge(defaultOptions, options);
        
        return hasRenderToArg ?
            new Highcharts.Chart(a, options, c) :
            new Highcharts.Chart(options, b);
    };

    const start = +new Date(),
    sparkline_tds = Array.from(document.querySelectorAll('td[data-sparkline]')),
    fullLen = sparkline_tds.length;
        
    let n = 0;
        
    // Creating 153 sparkline charts is quite fast in modern browsers, but IE8 and mobile
    // can take some seconds, so we split the input into chunks and apply them in timeouts
    // in order avoid locking up the browser process and allow interaction.
    function doChunk() {
        const time = +new Date(),
        len = sparkline_tds.length;
        
        for (let i = 0; i < len; i += 1) {
            const td = sparkline_tds[i];
            const stringdata = td.dataset.sparkline;
            const arr = stringdata.split('; ');
            const data = arr[0].split(', ').map(parseFloat);
            const chart = {};
        
            if (arr[1]) {
                chart.type = arr[1];
            }
        
            Highcharts.SparkLine(td, {
            series: [{
                data: data,
                pointStart: 1
            }],
            tooltip: {
                headerFormat: '<span style="font-size: 10px">' + td.parentElement.querySelector('th').innerText + ', Q{point.x}:</span><br/>',
                pointFormat: '<b>{point.y}.000</b> USD'
            },
            chart: chart
            });
        
            n += 1;

        }
    }
    doChunk();
}


$(document).ready(function () {
    var quant_table = document.getElementById("quantification_of_peptides");
    var thead = quant_table.getElementsByTagName("thead")[0];
    var tr = thead.getElementsByTagName("tr")[0];
    pep_header_ths = tr.getElementsByTagName("th");
    for(i=0; i<pep_header_ths.length; i++){
        $(pep_header_ths[i]).unbind("click");
    };
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

    draw_sparkline();

    $(".col-condition-sparkline").css("display", "none");

    $("#quantification_of_peptides .header").click(async function() {
        if(this.getAttribute("class").indexOf("rowheader") != -1){
            var span = "ID";
        } else{
            var span = this.getElementsByTagName("span")[0].innerText;
        };
        
        if((this.getAttribute("class").indexOf("headerSortUp") == -1) && (this.getAttribute("class").indexOf("headerSortDown") == -1)){
            $(this).attr('class', $(this).attr('class') + " headerSortUp");
            sortOrder = "headerSortUp";
        } else if(this.getAttribute("class").indexOf("headerSortUp") != -1) {
            $(this).attr('class', $(this).attr('class').replace("headerSortUp", "headerSortDown"));
            sortOrder = "headerSortDown";
        } else{
            $(this).attr('class', $(this).attr('class').replace("headerSortDown", "headerSortUp"));
            sortOrder = "headerSortUp";
        };
        
        for(i=0; i<pep_header_ths.length; i++){
            if(i==0){
                previous_span = "ID"
            } else {
                previous_span = pep_header_ths[i].getElementsByTagName("span")[0].innerText;
            }
            
            if(previous_span != span){
                pep_header_ths[i].setAttribute('class', pep_header_ths[i].className.replace(" headerSortUp", "").replace(" headerSortDown", ""));
            };
        };

        await quantFirst(sortOrder, span);
    });


    $("#peptide-distribution-button").click(function() {
        if(this.innerText == "distribution"){
            this.innerText = "intensity";        
            $("#quantification_of_peptides .col-condition").css("display", "none");
            $("#quantification_of_peptides .col-condition-sparkline").css("display", "table-cell");
        } else{
            this.innerText = "distribution";
            $("#quantification_of_peptides .col-condition").css("display", "table-cell");
            $("#quantification_of_peptides .col-condition-sparkline").css("display", "none");
        }

    });

})


//NextPage
async function quantNext(order, column){
    order = order||'original';
    currentRow = quantPageSize * quantPage;
    maxRow = currentRow + quantPageSize;
    if ( maxRow > numberRowsInQuantTable ) maxRow = numberRowsInQuantTable;
    await updateQuantData(currentRow, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = quant_trs[i].getElementsByTagName("td");
			quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
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
    
    quantPage++;

    if ( maxRow == numberRowsInQuantTable ) { quantNextText(); quantLastText(); }
    showPage(quantPageNum, quantPage);
    quantPreLink();
    quantFirstLink();
    draw_sparkline(); 
}

//PreviousPage
async function quantPre(order, column){
    order = order||'original';
    quantPage--;
    currentRow = quantPageSize * quantPage;
    maxRow = currentRow - quantPageSize;
    if ( currentRow > numberRowsInQuantTable ) currentRow = numberRowsInQuantTable;
	await updateQuantData(currentRow - quantPageSize, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = quant_trs[i].getElementsByTagName("td");
			quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}

        for (k = 1; k < 50; k++){
			quant_trs[k].style.display = '';
		}
	});

    if ( maxRow === 0 ){ quantPreText(); quantFirstText(); }
    showPage(quantPageNum, quantPage);
    quantNextLink();
    quantLastLink(); 
    draw_sparkline();
}

//FirstPage
async function quantFirst(order, column){
    order = order||'original';
    quantPage = 1;
    await updateQuantData(0, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = quant_trs[i].getElementsByTagName("td");
            quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
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
    draw_sparkline();
}

//LastPage
async function quantLast(order, column){
    order = order||'original';
    quantPage = parseInt(quantLastRows / quantPageSize + 1);
    await updateQuantData(quantLastRows, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = quant_trs[i].getElementsByTagName("td");
            quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
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
    draw_sparkline();
}

function showPage(pageNum, page){
    pageNum.innerHTML = page;

}

//TotalPage
async function quantPageCount(){
    let t;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
        .then(function (response) {
        let db = new window.SQL.Database(new Uint8Array(response.data));
        // execute query
        let s = new Date().getTime();
        let r = db.exec("select count(*) from PEPQUANT");
        let e = new Date().getTime();
        console.info("Time consuming to query data：" + (e - s) + "ms");
        t = r[0]['values'][0][0];
        })
        .catch(function (error) {
            console.info(error);
    });
    return t;
}

//TotalPage
async function updateQuantData(currentRow, order, column){
	let d;
    let r;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*',
			'Access-Control-Allow-Headers': '*'}})
                .then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));
			// query
			let s = new Date().getTime();
            if(order == "original"){
                r = db.exec("select * from PEPQUANT " + "limit "+ String(currentRow) + ",50" );
            } else if(order == "headerSortUp"){
                r = db.exec("select * from PEPQUANT " + "ORDER BY \"" + column + "\" DESC limit "+ String(currentRow) + ",50");
            } else{
                r = db.exec("select * from PEPQUANT " + "ORDER BY \"" + column + "\" ASC limit "+ String(currentRow) + ",50");
            }
			
			let e = new Date().getTime();
			console.info("Time consuming to query data：" + (e - s) + "ms");
			// parse data
			console.info(r);
			d = r[0]['values'];
    })
    .catch(function (error) {
        console.info(error);
        alert("Please set your browser to allow cross-domain requests.");
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

async function quant_page_jump(order, column){
	if(event.keyCode === 13){
        order = order||'original';
		quantPage = document.getElementById("pep_page").value;
		if (quantPage > parseInt(numberRowsInQuantTable / 50 + 1) || quantPage === ""){
			alert("not valid page!");
		} else if(quantPage === parseInt(numberRowsInQuantTable / 50 + 1)){
			await quantLast(order, column);
		} else{
			currentRow = quantPageSize * quantPage;
			maxRow = currentRow - quantPageSize;
			await updateQuantData(currentRow - quantPageSize, order, column).then(res =>{
                for(i=0; i<res.length; i++){
                    console.log(res[i]);
                    tds = quant_trs[i].getElementsByTagName("td");
                    quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
                    for(j=0; j<tds.length; j++){
                        if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                            tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                        } else if(res[i][j + 1] == null){
                            tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                        }
                        else{
                            tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                        }
                    }
                }
                console.log(i);
                for (k = 1; k < 50; k++){
                    if (k>=i) {
                        quant_trs[k].style.display = 'none';
                    } else{
                        quant_trs[k].style.display = '';
                    }
                }

            });

			if ( maxRow === 0 ){ quantPreText(); quantFirstText(); }
			showPage(quantPageNum, quantPage);
			quantPreLink();
			quantNextLink();
			quantFirstLink();
			quantLastLink();
            draw_sparkline();
		}
	}
}

async function searchData(filter, col, table){
	let d;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
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
        alert("Please set your browser to allow cross-domain requests.");
    });
	return d;
}


async function searchQuantFunction() {
    if (event.keyCode === 13) {
        var myInput=document.getElementById("quant_search");
        var filter=myInput.value.toUpperCase();
        var search_col=document.getElementById("quant_search_col");
        var index = search_col.selectedIndex;
        var value = search_col.options[index].text;

        await searchData(filter, value, 'PEPQUANT').then(res =>{
        for(i=0; i<res.length; i++){
            if(i>=50){
                break;
            }
            quant_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            tds = quant_trs[i].getElementsByTagName("td");
            for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
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
        draw_sparkline();
    }
}


// Protein Quantification Table
$(document).ready(function () {
    var prot_table = document.getElementById("quantification_of_protein");
    protTable = prot_table.getElementsByTagName("tbody")[0];
    var thead = prot_table.getElementsByTagName("thead")[0];
    var tr = thead.getElementsByTagName("tr")[0];
    header_ths = tr.getElementsByTagName("th");
    for(i=0; i<header_ths.length; i++){
        $(header_ths[i]).unbind("click");
    };

    protTotalPage = document.getElementById("protTotalPage");
    protPageNum = document.getElementById("protPageNum");

    protPreDom = document.getElementById("protPre");
    protNextDom = document.getElementById("protNext");
    protFirstDom = document.getElementById("protFirst");
    protLastDom = document.getElementById("protLast");
    prot_numrows = document.getElementById("quantification_of_protein_numrows_text");
    prot_sub = prot_numrows.getElementsByTagName("sub")[0];
    prot_trs = protTable.getElementsByTagName("tr");

    numberRowsInProtTable = 0;
    protPageSize = 50;
    protPage = 1;

    protPageCount().then(res => { protTotalPage.innerHTML = parseInt(res / 50 + 1);
        prot_sub.innerHTML = res;
        numberRowsInProtTable = res;
        protLastRows = protPageSize * (parseInt(res / 50 + 1) - 1);
    });
    protNextLink();
    protLastLink();
    protPageNum.innerHTML = '1';

    $("#quantification_of_protein .header").click(async function() {
        if(this.getAttribute("class").indexOf("rowheader") != -1){
            var span = "ProteinName";
        } else {
            var span = this.getElementsByTagName("span")[0].innerText;
        };
        
        if((this.getAttribute("class").indexOf("headerSortUp") == -1) && (this.getAttribute("class").indexOf("headerSortDown") == -1)){
            $(this).attr('class', $(this).attr('class') + " headerSortUp");
            sortOrder = "headerSortUp";
        } else if(this.getAttribute("class").indexOf("headerSortUp") != -1) {
            $(this).attr('class', $(this).attr('class').replace("headerSortUp", "headerSortDown"));
            sortOrder = "headerSortDown";
        } else{
            $(this).attr('class', $(this).attr('class').replace("headerSortDown", "headerSortUp"));
            sortOrder = "headerSortUp";
        };

        for(i=0; i<header_ths.length; i++){
            if(i==0){
                previous_span = "ProteinName";
            }else{
                previous_span = header_ths[i].getElementsByTagName("span")[0].innerText;
            };

            if(previous_span != span){
                header_ths[i].setAttribute('class', header_ths[i].className.replace(" headerSortUp", "").replace(" headerSortDown", ""));
            };
        };

        await protFirst(sortOrder, span);
    });
    
    $("#protein-distribution-button").click(function() {
        if(this.innerText == "distribution"){
            this.innerText = "intensity";        
            $("#quantification_of_protein .col-condition").css("display", "none");
            $("#quantification_of_protein .col-condition-sparkline").css("display", "table-cell");
        } else{
            this.innerText = "distribution";
            $("#quantification_of_protein .col-condition").css("display", "table-cell");
            $("#quantification_of_protein .col-condition-sparkline").css("display", "none");
        }

    });

})


//NextPage
async function protNext(order, column){
    order = order||'original';
    currentRow = protPageSize * protPage;
    maxRow = currentRow + protPageSize;
    if ( maxRow > numberRowsInProtTable ) maxRow = numberRowsInProtTable;
        await updateProtData(currentRow, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = prot_trs[i].getElementsByTagName("td");
			prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
        if ( i <= prot_trs.length - 1){
			for (k = i; k < prot_trs.length; k++){
				prot_trs[k].style.display = 'none';
			}
		}
	});
    protPage++;

    if ( maxRow == numberRowsInProtTable ) { protNextText(); protLastText(); }
    showPage(protPageNum, protPage);
    draw_sparkline();
    protPreLink();
    protFirstLink();
}

//PreviousPage
async function protPre(order, column){
    order = order||'original';
    protPage--;
    currentRow = protPageSize * protPage;
    maxRow = currentRow - protPageSize;
    if ( currentRow > numberRowsInProtTable ) currentRow = numberRowsInProtTable;
	await updateProtData(currentRow - protPageSize, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = prot_trs[i].getElementsByTagName("td");
			prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
			}
		}
        for (k = 1; k < 50; k++){
			prot_trs[k].style.display = '';
		}
	});


    if ( maxRow === 0 ){ protPreText(); protFirstText(); }
    showPage(protPageNum, protPage);
    draw_sparkline();
    protNextLink();
    protLastLink();
}

//FirstPage
async function protFirst(order, column){
    order = order||'original';
    protPage = 1;
    await updateProtData(0, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = prot_trs[i].getElementsByTagName("td");
            prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
            }
        }
		console.log(i);
		for (k = 1; k < 50; k++){
			prot_trs[k].style.display = '';
		}

    });
    showPage(protPageNum, protPage);
    draw_sparkline();
    protPreText();
    protNextLink();
    protLastLink();
}

//LastPage
async function protLast(order, column){
    order = order||'original';
    protPage = parseInt(protLastRows / protPageSize + 1);
    await updateProtData(protLastRows, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = prot_trs[i].getElementsByTagName("td");
            prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                    tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                } else if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				}
				else{
					tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
				}
            }
        }
		if ( i <= prot_trs.length - 1){
			for (k = i; k < prot_trs.length; k++){
                prot_trs[k].style.display = 'none';
			}
		}
    });
    showPage(protPageNum, protPage);
    draw_sparkline();
    protPreLink();
    protNextText();
    protFirstLink();
}

//TotalPage
async function protPageCount(){
    let t;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
        .then(function (response) {
        let db = new window.SQL.Database(new Uint8Array(response.data));
        let r = db.exec("select count(*) from PROTQUANT");
        t = r[0]['values'][0][0];
        })
        .catch(function (error) {
            console.info(error);
    });
    return t;
}

//ShowLink
function protPreLink(){ protPreDom.innerHTML = "<a href='javascript:protPre();'>Previous Page</a>";}
function protPreText(){ protPreDom.innerHTML = "Previous Page";}

function protNextLink(){ protNextDom.innerHTML = "<a href='javascript:protNext();'>Next Page</a>";}
function protNextText(){ protNextDom.innerHTML = "Next Page";}

function protFirstLink(){ protFirstDom.innerHTML = "<a href='javascript:protFirst();'>First Page</a>";}
function protFirstText(){ protFirstDom.innerHTML = "First Page";}

function protLastLink(){ protLastDom.innerHTML = "<a href='javascript:protLast();'>Last Page</a>";}
function protLastText(){ protLastDom.innerHTML = "Last Page";}

async function updateProtData(currentRow, order, column){
	let d;
    let r;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
    		.then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));
            if(order == "original"){
                r = db.exec("select * from PROTQUANT " + "limit "+ String(currentRow) + ",50");
            } else if(order == "headerSortUp"){
                r = db.exec("select * from PROTQUANT " + "ORDER BY \"" + column + "\" DESC limit "+ String(currentRow) + ",50");
            } else{
                r = db.exec("select * from PROTQUANT " + "ORDER BY \"" + column + "\" ASC limit "+ String(currentRow) + ",50");
            }
			d = r[0]['values'];
    })
    .catch(function (error) {
        console.info(error);
        alert("Please set your browser to allow cross-domain requests.");
    });
	return d;
}

async function prot_page_jump(order, column){
	if(event.keyCode === 13){
        order = order||'original';
		protPage = document.getElementById("prot_page").value;
		if (protPage > parseInt(numberRowsInProtTable / 50 + 1) || protPage === ""){
			alert("not valid page!");
		}

		else if(protPage === parseInt(numberRowsInProtTable / 50 + 1)){
			await protLast(order, column);
		} else{
			currentRow = protPageSize * protPage;
			maxRow = currentRow - protPageSize;
			await updateProtData(currentRow - protPageSize, order, column).then(res =>{
                for(i=0; i<res.length; i++){
                    tds = prot_trs[i].getElementsByTagName("td");
                    prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
                    for(j=0; j<tds.length; j++){
                        if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                            tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                        } else if(res[i][j + 1] == null){
                            tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                        }
                        else{
                            tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                        }
                    }
                }
                for (k = 1; k < 50; k++){
                    if (k>=i) {
                        prot_trs[k].style.display = 'none';
                    } else{
                        prot_trs[k].style.display = '';
                    }
                }
            });

			if ( maxRow === 0 ){ protPreText(); protFirstText(); }
			showPage(protPageNum, protPage);
            draw_sparkline();
			protPreLink();
			protNextLink();
			protFirstLink();
			protLastLink();
		}
	}
}

async function searchProtFunction() {
    if (event.keyCode === 13) {
        var myInput=document.getElementById("prot_search");
        var filter=myInput.value.toUpperCase();
        var search_col=document.getElementById("prot_search_col");
        var index = search_col.selectedIndex;
        var value = search_col.options[index].text;

        await searchData(filter, value, 'PROTQUANT').then(res =>{
            for(i=0; i<res.length; i++){
                if(i>=50){
                    break;
                }
                prot_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
                tds = prot_trs[i].getElementsByTagName("td");
                for(j=0; j<tds.length; j++){
                    if(tds[j].getAttribute("class") == "data-sparkline col-condition-sparkline"){
                        tds[j].setAttribute("data-sparkline", String(res[i][j + 1]));
                    } else if(res[i][j + 1] == null){
                        tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                    }
                    else{
                        tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                    }
                }

            }
            for (k = 0; k < 50; k++){
                if (k>=i) {
                    prot_trs[k].style.display = 'none';
                } else{
                    prot_trs[k].style.display = '';
                }
            }
        });
        draw_sparkline();
    }
}