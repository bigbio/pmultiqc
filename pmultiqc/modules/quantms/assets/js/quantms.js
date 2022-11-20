function draw_sparkline(table_dict){
    Highcharts.SparkLine = function (a, b, c) {
        const hasRenderToArg = typeof a === 'string' || a.nodeName;
        let options = arguments[hasRenderToArg ? 1 : 0];
        const defaultOptions = {
            chart: {
                renderTo: (options.chart && options.chart.renderTo) || (hasRenderToArg && a),
                backgroundColor: null,
                borderWidth: 0,
                type: 'area',
                margin: [0, 0, 0, 0],
                //width: 120,
                //height: 40,
                style: {
                    overflow: 'visible'
                },
                // small optimization, saves 1-2 ms each sparkline
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
                tickPositions: [],
                lineWidth: 1,
                lineColor: '#C0D0E0',
                gridLineWidth: 0
            },
            yAxis: {
                min: 0.0,
                //max: table_dict["maxValue"],
                endOnTick: false,
                startOnTick: false,
                labels: {
                    enabled: false
                },
                title: {
                    text: null
                },
                tickPositions: [0],
                gridLineWidth: 0
            },
            legend: {
                enabled: false
            },
            dataLabels: {
                enabled: true,
                color: '#000000',
                nullFormat: 'N/A'
            },
            tooltip: {
                hideDelay: 0,
                outside: true,
                nullFormat: 'N/A',
                backgroundColor: '#fff'
            },
            exporting: {
                enabled: false
            },
            plotOptions: {
                series: {
                    enableMouseTracking: true,
                    animation: false,
                    lineWidth: 1,
                    minPointLength: 2,
                    shadow: false,
                    states: {
                        hover: {
                            lineWidth: 3
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

    sparkline_tds = Array.from(document.querySelectorAll(table_dict['name'] + ' ' + 'td[data-sparkline]'));
    sparkline_tds_cell = Array.from(document.querySelectorAll(table_dict['name'] + ' ' + 'td.data-sparkline'));
    len_tds_tr = Array.from(document.querySelectorAll('th.col-condition-sparkline')).length;
    average_intensity_col = Array.from(document.querySelectorAll(table_dict['name'] + ' ' + 'td.Average_Intensity > .wrapper > .val'));
    if(average_intensity_col.length == 0){
        average_intensity_col = Array.from(document.querySelectorAll(table_dict['name'] + ' ' + 'td.Average_Spectrum_Counting > .wrapper > .val'));
    }

    function doChunk() {
        len = sparkline_tds.length;
        
        for (let i = 0; i < len; i += 1) {
            const td = sparkline_tds[i];
            const stringdata = td.dataset.sparkline;
            //TODO figure out when None and when nan happens
            if (stringdata === "nan" || stringdata === "None") continue;
            // decode commas etc
            const data = JSON.parse(decodeURIComponent(stringdata))
            const series = []
            for (const sample in data)
            {
                series.push({data: [data[sample]], name: sample})
            }
            const chart = {};
            chart.type = 'column'
        
            td.setAttribute("width", series.length * 30)
            chart.width = td.getAttribute('width')
            chart.height = 50

            Highcharts.SparkLine(td, {
            series: series,
            tooltip: {
                padding: 2,
                style:
                {
                    fontSize: "10px"
                },
                outside: true,
                headerFormat: '{series.name}: ',
                pointFormat: '<b>{point.y}</b>'
            },
            chart: chart
            });

            average_intensity = parseFloat(average_intensity_col[parseInt(i / len_tds_tr)].innerText);
            rects = sparkline_tds_cell[i].querySelectorAll("svg > .highcharts-series-group > .highcharts-series > rect");

            for(let g = 0; g < rects.length; g++){
                if(parseFloat(data[g]) < average_intensity){
                    rects[g].setAttribute("fill", "#910000");
                }
            };
        }
    }
    doChunk();
}


$(document).ready(function () {
    $(".col-condition-sparkline").css("display", "none");
    var quant_table = document.getElementById("quantification_of_peptides");
    if (quant_table == null){
        return
    }
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
    quantPageCount().then(res => { quantTotalPage.innerHTML = parseInt(res / 51 + 1);
        quant_sub.innerHTML = res;
        numberRowsInQuantTable = res;
        quantLastRows = quantPageSize * (parseInt(res / 50 + 1) - 1);
    });
    quantNextLink();
    quantLastLink();
    quantPageNum.innerHTML = '1';
    pep_maxValue = quant_table.getAttribute("data-dmax");
    peptide_table_dict = {'name': '#quantification_of_peptides', 'maxValue': pep_maxValue};

    draw_sparkline(peptide_table_dict);

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
        if(this.innerText == " Show replicates"){
            $("#quantification_of_peptides tr").css("height", "50px");
            this.innerHTML = "<span class='glyphicon glyphicon glyphicon-stats'></span> Hide replicates";        
            $("#quantification_of_peptides .col-condition").css("display", "none");
            $("#quantification_of_peptides .col-condition-sparkline").css("display", "table-cell");
        } else{
            $("#quantification_of_peptides tr").css("height", "100%");
            this.innerHTML = "<span class='glyphicon glyphicon glyphicon-stats'></span> Show replicates"; 
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
				} else if(j == 3){
					tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
    draw_sparkline(peptide_table_dict); 
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
				} else if(j == 3){
					tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
    draw_sparkline(peptide_table_dict);
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
				} else if(j == 3){
					tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
    draw_sparkline(peptide_table_dict);
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
				} else if(j == 3){
					tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
    draw_sparkline(peptide_table_dict);
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
                        } else if(j == 3){
                            tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
            draw_sparkline(peptide_table_dict);
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
				} else if(j == 3){
					tds[j].getElementsByClassName('val')[0].innerHTML = exponential_form(res[i][j + 1]);
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
        draw_sparkline(peptide_table_dict);
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
    prot_maxValue = prot_table.getAttribute("data-dmax");
    protein_table_dict = {'name': '#quantification_of_protein', 'maxValue': prot_maxValue};

    draw_sparkline(protein_table_dict);

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
        if(this.innerText == " Show replicates"){
            $("#quantification_of_protein tr").css("height", "50px");
            this.innerHTML = "<span class='glyphicon glyphicon glyphicon-stats'></span> Hide replicates";        
            $("#quantification_of_protein .col-condition").css("display", "none");
            $("#quantification_of_protein .col-condition-sparkline").css("display", "table-cell");
        } else{
            $("#quantification_of_protein tr").css("height", "100%");
            this.innerHTML = "<span class='glyphicon glyphicon glyphicon-stats'></span> Show replicates"; 
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
    draw_sparkline(protein_table_dict);
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
    draw_sparkline(protein_table_dict);
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
		for (k = 1; k < i; k++){
			prot_trs[k].style.display = '';
		}

    });
    showPage(protPageNum, protPage);
    draw_sparkline(protein_table_dict);
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
    draw_sparkline(protein_table_dict);
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
            draw_sparkline(protein_table_dict);
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
        draw_sparkline(protein_table_dict);
    }
}


// PSM Table
$(document).ready(function () {
    var psm_table = document.getElementById("peptide_spectrum_matches");
    if (psm_table == null){
        return
    }
    psmTable = psm_table.getElementsByTagName("tbody")[0];
    var thead = psm_table.getElementsByTagName("thead")[0];
    var tr = thead.getElementsByTagName("tr")[0];
    header_ths = tr.getElementsByTagName("th");
    for(i=0; i<header_ths.length; i++){
        $(header_ths[i]).unbind("click");
    };

    psmTotalPage = document.getElementById("psmTotalPage");
    psmPageNum = document.getElementById("psmPageNum");

    psmPreDom = document.getElementById("psmPre");
    psmNextDom = document.getElementById("psmNext");
    psmFirstDom = document.getElementById("psmFirst");
    psmLastDom = document.getElementById("psmLast");
    psm_numrows = document.getElementById("peptide_spectrum_matches_numrows_text");
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
    psm_maxValue = psm_table.getAttribute("data-dmax");

    $("#peptide_spectrum_matches .header").click(async function() {
        if(this.getAttribute("class").indexOf("rowheader") != -1){
            var span = "PSM_ID";
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
                previous_span = "PSM_ID";
            }else{
                previous_span = header_ths[i].getElementsByTagName("span")[0].innerText;
            };

            if(previous_span != span){
                header_ths[i].setAttribute('class', header_ths[i].className.replace(" headerSortUp", "").replace(" headerSortDown", ""));
            };
        };

        await psmFirst(sortOrder, span);
    });

})


//NextPage
async function psmNext(order, column){
    order = order||'original';
    currentRow = psmPageSize * psmPage;
    maxRow = currentRow + psmPageSize;
    if ( maxRow > numberRowsInPsmTable ) maxRow = numberRowsInPsmTable;
        await updatePsmData(currentRow, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = psm_trs[i].getElementsByTagName("td");
			psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				} else if (j==3){
					tds[j].innerHTML = exponential_form(res[i][j + 1]);
				} else{
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
    psmPage++;

    if ( maxRow == numberRowsInPsmTable ) { psmNextText(); psmLastText(); }
    showPage(psmPageNum, psmPage);
    psmPreLink();
    psmFirstLink();
}

//PreviousPage
async function psmPre(order, column){
    order = order||'original';
    psmPage--;
    currentRow = psmPageSize * psmPage;
    maxRow = currentRow - psmPageSize;
    if ( currentRow > numberRowsInPsmTable ) currentRow = numberRowsInPsmTable;
	await updatePsmData(currentRow - psmPageSize, order, column).then(res =>{
		for(i=0; i<res.length; i++){
			console.log(res[i]);
			tds = psm_trs[i].getElementsByTagName("td");
			psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
			for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				} else if (j==3){
					tds[j].innerHTML = exponential_form(res[i][j + 1]);
				} else{
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
    psmNextLink();
    psmLastLink();
}

//FirstPage
async function psmFirst(order, column){
    order = order||'original';
    psmPage = 1;
    await updatePsmData(0, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = psm_trs[i].getElementsByTagName("td");
            psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				} else if (j==3){
					tds[j].innerHTML = exponential_form(res[i][j + 1]);
				} else{
                    tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                }
            }
        }
		console.log(i);
		for (k = 1; k < i; k++){
			psm_trs[k].style.display = '';
		}

    });
    showPage(psmPageNum, psmPage);
    psmPreText();
    psmNextLink();
    psmLastLink();
}

//LastPage
async function psmLast(order, column){
    order = order||'original';
    psmPage = parseInt(psmLastRows / psmPageSize + 1);
    await updatePsmData(psmLastRows, order, column).then(res =>{
        for(i=0; i<res.length; i++){
            console.log(res[i]);
            tds = psm_trs[i].getElementsByTagName("td");
            psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
            for(j=0; j<tds.length; j++){
                if(res[i][j + 1] == null){
					tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
				} else if (j==3){
					tds[j].innerHTML = exponential_form(res[i][j + 1]);
				} else{
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
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
        .then(function (response) {
        let db = new window.SQL.Database(new Uint8Array(response.data));
        let r = db.exec("select count(*) from PSM");
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

async function updatePsmData(currentRow, order, column){
	let d;
    let r;
    await axios.get("quantms.db", {responseType: 'arraybuffer'}, {headers:{'Access-Control-Allow-Origin': '*'}})
    		.then(function (response) {
			let db = new window.SQL.Database(new Uint8Array(response.data));
            if(order == "original"){
                r = db.exec("select * from PSM " + "limit "+ String(currentRow) + ",50");
            } else if(order == "headerSortUp"){
                r = db.exec("select * from PSM " + "ORDER BY \"" + column + "\" DESC limit "+ String(currentRow) + ",50");
            } else{
                r = db.exec("select * from PSM " + "ORDER BY \"" + column + "\" ASC limit "+ String(currentRow) + ",50");
            }
			d = r[0]['values'];
    })
    .catch(function (error) {
        console.info(error);
        alert("Please set your browser to allow cross-domain requests.");
    });
	return d;
}

async function psm_page_jump(order, column){
	if(event.keyCode === 13){
        order = order||'original';
		psmPage = document.getElementById("psm_page").value;
		if (psmPage > parseInt(numberRowsInPsmTable / 50 + 1) || psmPage === ""){
			alert("not valid page!");
		}

		else if(psmPage === parseInt(numberRowsInPsmTable / 50 + 1)){
			await psmLast(order, column);
		} else{
			currentRow = psmPageSize * psmPage;
			maxRow = currentRow - psmPageSize;
			await updatePsmData(currentRow - psmPageSize, order, column).then(res =>{
                for(i=0; i<res.length; i++){
                    tds = psm_trs[i].getElementsByTagName("td");
                    psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
                    for(j=0; j<tds.length; j++){
                        if(res[i][j + 1] == null){
                            tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                        } else if (j==3){
                            tds[j].innerHTML = exponential_form(res[i][j + 1]);
                        } else{
                            tds[j].getElementsByClassName('val')[0].innerHTML = res[i][j + 1];
                        }
                    }
                }
                for (k = 1; k < 50; k++){
                    if (k>=i) {
                        psm_trs[k].style.display = 'none';
                    } else{
                        psm_trs[k].style.display = '';
                    }
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

async function searchPsmFunction() {
    if (event.keyCode === 13) {
        var myInput=document.getElementById("psm_search");
        var filter=myInput.value.toUpperCase();
        var search_col=document.getElementById("psm_search_col");
        var index = search_col.selectedIndex;
        var value = search_col.options[index].text;

        await searchData(filter, value, 'PSM').then(res =>{
            for(i=0; i<res.length; i++){
                if(i>=50){
                    break;
                }
                psm_trs[i].getElementsByTagName("th")[0].innerHTML = res[i][0];
                tds = psm_trs[i].getElementsByTagName("td");
                for(j=0; j<tds.length; j++){
                    if(res[i][j + 1] == null){
                        tds[j].getElementsByClassName('val')[0].innerHTML = String(res[i][j + 1]);
                    } else if (j==3){
                        tds[j].innerHTML = exponential_form(res[i][j + 1]);
                    } else{
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

$(document).ready(function () {
    var search_1 = $('[data-target="search_scores_summary-1"]').parent();
    var search_2 = $('[data-target="search_scores_summary-2-1"]').parent();
    var pep = $('[data-target="search_engine_PEP-1"]').parent();
    var consensus = $('[data-target="consensus_summary-1"]').parent();
    if(search_1.length > 0){
        search_1.prepend("<span id='search_btn' class='btn btn-default btn-sm'>select files</span>");
        search_1.children().css({"position":"inherit","float":"inherit","width":"200px","overflow-x":"auto"});
        search_1.css({"height": "30px","overflow-y":"hidden"});
        search_1.mouseenter(function () {
            $(this).children("button").show();
            $(this).css({"height": "140px","overflow-y":"auto"});
        });
        search_1.mouseleave(function () {
            $(this).children("button").hide();
            $(this).css({"height": "30px","overflow-y":"auto"});
        });
        search_1.children().each(function () {
            $(this).css({"display": "block"});
        });
    };
    if(search_2.length > 0){
        search_2.prepend("<span id='search_btn' class='btn btn-default btn-sm'>select files</span>");
        search_2.children().css({"position":"inherit","float":"inherit","width":"200px","overflow-x":"auto"});
        search_2.css({"height": "30px","overflow-y":"hidden"});
        search_2.mouseenter(function () {
            $(this).children("button").show();
            $(this).css({"height": "140px","overflow-y":"auto"});
        });
        search_2.mouseleave(function () {
            $(this).children("button").hide();
            $(this).css({"height": "30px","overflow-y":"auto"});
        });
        search_2.children().each(function () {
            $(this).css({"display": "block"});
        });
    };
    if(pep.length > 0){
        pep.prepend("<span id='search_btn' class='btn btn-default btn-sm'>select files</span>");
        pep.children().css({"position":"inherit","float":"inherit","width":"200px","overflow-x":"auto"});
        pep.css({"height": "30px","overflow-y":"hidden"});
        pep.mouseenter(function () {
            $(this).children("button").show();
            $(this).css({"height": "140px","overflow-y":"auto"});
        });
        pep.mouseleave(function () {
            $(this).children("button").hide();
            $(this).css({"height": "30px","overflow-y":"auto"});
        });
        pep.children().each(function () {
            $(this).css({"display": "block"});
        });
    };
    if(consensus.length > 0){
        consensus.prepend("<span id='search_btn' class='btn btn-default btn-sm'>select files</span>");
        consensus.children().css({"position":"inherit","float":"inherit","width":"200px","overflow-x":"auto"});
        consensus.css({"height": "30px","overflow-y":"hidden"});
        consensus.mouseenter(function () {
            $(this).children("button").show();
            $(this).css({"height": "140px","overflow-y":"auto"});
        });
        consensus.mouseleave(function () {
            $(this).children("button").hide();
            $(this).css({"height": "30px","overflow-y":"auto"});
        });
        consensus.children().each(function () {
            $(this).css({"display": "block"});
        });
    };
})

function exponential_form(value){
    if (value == null || value == '0.00e+00' || isNaN(value)){
        return String(value)
    } else{
        var p = Math.floor(Math.log(value)/Math.LN10);
        var n = value * Math.pow(10, -p);
        n = n.toFixed(5)
        return n + 'e' + p
    }
}
