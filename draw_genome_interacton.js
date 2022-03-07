function get_locus_xy(locus, chr_len, chr_x, chr_y, chr_radius) {
    var theta = locus * 2 * Math.PI / chr_len;
    
    var xx = chr_radius * Math.sin(theta) + chr_x;

    var yy = chr_radius * Math.cos(theta) + chr_y;

    var new_theta = 360 * locus / chr_len;

  return [xx, yy, new_theta]
}

function add_arrow(s, tag){
    let draw_loci_fig = tag['fig'];
    let x1 = draw_loci_fig.x1;
    let y1 = draw_loci_fig.y1;
    let new_theta = draw_loci_fig.new_theta;
    let strand = tag.loci_strand;
    new_theta = 90-new_theta;
    let arrow_color = "#b3461e";
    if(strand == "-"){
        new_theta += 180;
        arrow_color = "blue";
    }
    let w = 6 //The width lower bound of arrow 
    let wa = 16//The weidth of arrow 
    let h = 25 // total height
    let ha = 6// arrow of height
    let jiancha = (wa-w)/2//The higher and lower place of wideth difference
    let hs = h - ha
    let loci = s.polyline([x1, y1, x1+w, y1, x1+w, y1-hs, x1+w+jiancha, y1-hs, x1+(w/2), y1-h, x1-jiancha, y1-hs, x1, y1-hs]);
    loci.attr({
        fill:arrow_color,
        stroke:"black",
        })
    let m = new Snap.Matrix();
    m.rotate(new_theta, x1,y1);
    loci.transform(m)
}

function add_text(s, x1, y1, tag, new_theta){
    let refseq = tag.nc_id;
    console.log("refseq");
    console.log(tag);
    let region_start = parseInt(tag.loci_start) - contigs[refseq].start;
    let region_end = parseInt(tag.loci_end) - contigs[refseq].start;
    let strand = tag.loci_strand;
    let x_offset = 20;
    let y_offset = 6;
    // if(new_theta>180){
    //     x_offset = -140;
    // }
    let text = s.text(x1+x_offset,y1+y_offset,"Acr-Aca "+tag.loci_id+": "+region_start+"-"+region_end+" ("+strand+")");
    let color = "blue";
    if (strand=="+"){
        color = "#b3461e";
    }
    text.attr({
            fill:color,
            fontSize:10
    })
}

function draw_loci(s, tag){
    let draw_loci_fig = tag['fig'];
    add_arrow(s, tag);
    add_text(s, draw_loci_fig.text_x1, draw_loci_fig.text_y1, tag, draw_loci_fig.new_theta);
}
function parse_5kbp(_5kbp){
    var mydict = {};
    var arr_5k = _5kbp.split(",");
    for(let i5=0; i5<arr_5k.length; i5++){
        mydict[i5] = {}
        let tmp_arr = arr_5k[i5].split("|");
        for(j5=0; j5<tmp_arr.length; j5++){
            if (tmp_arr[j5].startsWith("CAS-Type") || tmp_arr[j5].startsWith("CAS(") ){
                cas_ids = tmp_arr[j5].split("+");
                let cases = {};
                for(z5=0; z5<cas_ids.length; z5++){
                    //CAS-TypeIC(457825-465657)
                    [id, loc] = cas_ids[z5].split("(");
                    [start, end] = loc.split(")")[0].split("-");
                    cases[z5] =  {'id':id, 'start':start, 'end':end};
                }
                mydict[i5]['cas'] = cases;
            }
            else{
                [key, value] = tmp_arr[j5].split("=");
                mydict[i5][key] = value;
            }
            
        }
    }
    return mydict;   
}


function add_cas(s, tt, centerX, centerY, chr_len, chromRadius){
    
    let cas = tt.sort(compare("start"));
    let text_y = 100000;
    let preTheta = 0;
    for(i in cas){
        one = cas[i];
        castype = one.id;
        cas_start = parseFloat(one.start);
        cas_end = parseFloat(one.end);
        cas_strand = one.strand;
        [x1, y1, new_theta]= get_locus_xy(cas_start, chr_len, centerX, centerY, chromRadius);
        //add rect
        var rec_cas = s.rect(x1-5,y1,5,5);
        rec_cas.attr({
            fill:"#d11ba1",
            stroke:"#d11ba1",
            strokeWidth:3
        }); 
        //add text
        if(new_theta<180){
            x_offset = -150;
        }else{
            x_offset =20;
        }
        cur_text_y = get_cur_texty(text_y, y1, preTheta,new_theta,20);
        var text2 = s.text(x1+x_offset,cur_text_y,castype+":"+cas_start+"-"+cas_end);
        
        text2.attr({
            fill:"#d11ba1",
            fontSize:10
        });
        text_y = cur_text_y;
    }
    preTheta = new_theta;

}

function get_cur_texty(preTextY, curTextY, preTheta, curTheta, gap){
    if(preTextY==100000){
        curTextY = curTextY;
    }
    else if(preTheta<180 && curTheta<180){
        if(curTextY < preTextY){
            if(preTextY  - curTextY > gap){
                curTextY = curTextY
            }
            else{
                curTextY -= gap;
            }
        }
        else{
               curTextY = preTextY - gap;
        }
    }
    else if(preTheta>180 && curTheta>180){
        if(curTextY > preTextY){
            if(curTextY - preTextY > gap){
                curTextY = curTextY
            }
            else{
                curTextY += gap;
            }
        }
        else{
               curTextY = preTextY + gap;
        }
    }
    return curTextY;
}
function add_shape(s, curElement, preElement, eleConfig, centerX, centerY, elementNumber, chr_len, chromRadius){
    // shape of the tag , strokeColor, fontColor, spreadDistance, tagName, judgeDistance, gap
    /*
    example for target
    eleConfig={
        "shape":"ellipse",
        "strokeColor":"green",
        "fontColor":"black",
        "spreadDistance":200,
        "tagName":"Target",
        "judgeDistance":2500,
        "gap": 20
    }
    */
    // eleConfig is obj  {}
    var curX, curY, curTheta;
    var curStart = curElement['curStart'], curEnd = curElement['curEnd'];
    var preTextX = preTextY = "100000", preTheta = "1";
    
    if(!jQuery.isEmptyObject(preElement)){
        preTextX = preElement['preTextX'];
        preTextY = preElement['preTextY'];
        preTheta = preElement['preTheta'];
    }

    [curX, curY, curTheta]= get_locus_xy(curStart, chr_len, centerX, centerY, chromRadius);
    
    //add shape
    if(eleConfig.shape == "ellipse"){
        var shape_fig = s.ellipse(curX, curY, 5, 10);
        shape_fig.attr({
            fill:"white",
            stroke:eleConfig.strokeColor,
        });
    }
    
    //add text
    // hint: spreadDistance: Target 200, spacer:100, acr-aca:50 
    var bigR = chromRadius + eleConfig.spreadDistance;
    var t_offset = 0;
    if(curTheta > 180){
        t_offset = -150;
    }
    var arc = Math.PI * curTheta / 180 //转成弧度制
    var curTextX = bigR * Math.sin(arc) + centerX;
    var curTextY = bigR * Math.cos(arc) + centerY;
    //dynamic addjust the text location
    
    curTextY = get_cur_texty(preTextY, curTextY, preTheta, curTheta, eleConfig.gap);
    
    //add target tag
    var textFig = s.text(curTextX + t_offset, curTextY, elementNumber.toString() + " "+eleConfig.tagName+" Position"+":"+curStart+"-"+curEnd);
    textFig.attr({
            fill: eleConfig.fontColor,
            fontSize:10
    });
    // draw line for connect target round and tag
    var lineFig = s.line(curX, curY, curTextX, curTextY);
    lineFig.attr({
        fill: "9CA4A2", //color to fill
        stroke: "9CA4A2", // color for the line
        strokeWidth: 1,//width of the line
        fillOpacity: 0.5 //Opacity
    });
    var pe = {
        "preTextX":curTextX,
        "preTextY":curTextY,
        "preTheta":curTheta
    }
    return pe;
}

function draw_5kbp_target(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius){
    var last_target_pos = "";
    var target_obj=[];
    var target_i=0;
    for(qwe in _5kbp_dict){
        temp = _5kbp_dict[qwe];
        
        if(temp.Target_Pos!=last_target_pos){
            var target_arr_temp = temp.Target_Pos.split("-");
            console.log(temp);
            target_obj[target_i] = {
                "curStart":target_arr_temp[0],
                "curEnd":target_arr_temp[1]
            };
            last_target_pos = temp.Target_Pos;
            target_i += 1;
        } 
    }
    var sort_Obj = target_obj.sort(compare("curStart"));
    var tti;
    var preElement = {};
    var eleConfig={
        "shape":"ellipse",
        "strokeColor":"#999900",
        "fontColor":"#999900",
        "spreadDistance":300,
        "tagName":"Target",
        "judgeDistance":2500,
        "gap": 20
    };
    for(tti=0; tti<sort_Obj.length; tti++){
        preElement = add_shape(s, sort_Obj[tti], preElement, eleConfig, centerX, centerY, tti+1, chr_len, chromRadius);
    }
}

function draw_5kbp_spacer(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius){
    var last_spacer_pos = "";
    var spacer_obj=[];
    var spacer_i=0;
    for(qwe in _5kbp_dict){
        temp = _5kbp_dict[qwe];
        
        if(temp.Spacer_Pos!=last_spacer_pos){
            var arr_temp = temp.Spacer_Pos.split("-");
            spacer_obj[spacer_i] = {
                "curStart":arr_temp[0],
                "curEnd":arr_temp[1]
            };
            last_spacer_pos = temp.Spacer_Pos;
            spacer_i += 1;
        }  
    }
    var sort_Obj = spacer_obj.sort(compare("curStart"));
    var tti;
    var preElement = {};
    var eleConfig={
        "shape":"ellipse",
        "strokeColor":"purple",
        "fontColor":"purple",
        "spreadDistance":150,
        "tagName":"Spacer",
        "judgeDistance":16000,
        "gap": 15
    };
    for(tti=0; tti<sort_Obj.length; tti++){
        preElement = add_shape(s, sort_Obj[tti], preElement, eleConfig, centerX, centerY, tti+1, chr_len, chromRadius)
    }
}

function draw_5kbp_target_spacer_connection(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius){

    var c = _5kbp_dict;
    for(var ci in c){
        var start_pos = c[ci].Spacer_Pos.split("-")[0];
        var end_pos = c[ci].Target_Pos.split("-")[0];
        [sxx, syy, theta1] = get_locus_xy(start_pos, chr_len, centerX, centerY, chromRadius);
        [exx, eyy, theta2] = get_locus_xy(end_pos, chr_len, centerX, centerY, chromRadius);
        sweep_flag="0";
        bigArcFlag = "0";
        if(parseInt(start_pos) < parseInt(end_pos)){
            if (Math.abs(theta2-theta1)>180){
                sweep_flag = "0";
            }
            else{
                sweep_flag="1";
            }
            
        }
        // if(Math.sqrt(Math.pow(sxx - exx, 2) + Math.pow(syy - eyy, 2)) > chromRadius){
        //     bigArcFlag = "1";
        // }   l; 
        var chord_len = Math.sqrt(Math.pow(sxx - exx, 2) + Math.pow(syy - eyy, 2));
        var arc_theta = Math.asin(chord_len/(2*chromRadius));
        var arc_radius = chord_len/(2*Math.sin(Math.PI/2 - arc_theta));
        let mypath = `M${sxx} ${syy} A${arc_radius} ${arc_radius}, 0, ${bigArcFlag} ${sweep_flag}, ${exx} ${eyy}`;
        // mypath = 
        // "M" + sxx + " " + syy+ " " +
        // "A" + arc_radius + " " + arc_radius +"," + " " + // arc_radius height and width
        // "0," + " " + bigArcFlag+" "+sweep_flag+"," + " " +  // theta, big arc flag, sweep flag
        // exx + " " + eyy;
        var arc_fig  =  s.path(mypath);
        arc_fig.attr({
            stroke:"blue",
            fill:"none"
        });

    }
}

function compare(property){
    return function(obj1, obj2){
        var value1 = obj1[property];
        var value2 = obj2[property];
        return value1 - value2;
    }
}

function add_strand_crispr(s, temppp, result_json, centerX, centerY, chr_len, chromRadius, crispr_cas_start_y){
    var _5kbp_dict = temppp;
    // add strand to cas 
    for(let i in _5kbp_dict){
        // look for which crispr cas pair
        for(let h in result_json.Sequences){
            let Version = result_json.Sequences[h].Version;
            if(_5kbp_dict[i]['Spacer Accession']===Version){
                console.log(_5kbp_dict[i]);
                let kk = result_json.Sequences[h];
                let cas_arr = kk.Cas;
                let crispr_arr = kk.Crisprs;
                for(let j in _5kbp_dict[i].cas){
                    //1. add strand to cas
                    
                    var c1 = _5kbp_dict[i].cas[j]; // c1 is 5kbp
                    for(let z in cas_arr){
                        var c2 = cas_arr[z]; // c2 is from json
                        if (c1.start == c2.Start && c1.id == c2.Type){
                            _5kbp_dict[i].cas[j]['strand']  = c2.Genes[0].Orientation;
                            _5kbp_dict[i].cas[j]['Genes']  = c2.Genes;
                            // console.log(c2);
                            break;
                        }
                    }
                    //2. add crispr
                    var spacer_pos = _5kbp_dict[i].Spacer_Pos;
                    [sstart, send] = spacer_pos.split("-");
                    for(let z in crispr_arr){
                        var crispr_one = crispr_arr[z];
                        if(crispr_one.Start<= sstart && crispr_one.End >= send){
                            _5kbp_dict[i]['crispr'] = crispr_one;
                            // console.log(crispr_one);
                            break;
                        }
                    }
                }
                break;
            }
            
            
            
        }
        
    }
    console.log("is assemble???");
    console.log(is_assemble);
    if(is_assemble==="1"){
        draw_5kbp_target(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius);
        draw_5kbp_spacer(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius);
        draw_5kbp_target_spacer_connection(s, _5kbp_dict, centerX, centerY, chr_len, chromRadius);
    }
    
    // add_cas(_5kbp_dict[0].cas, centerX, centerY);
    s.text(20,crispr_cas_start_y, "CRISPR-Cas Array (Self Targeting Spacers)");
    console.log(_5kbp_dict);
    add_crispr_cas(s, _5kbp_dict, result_json, crispr_cas_start_y, centerX, centerY, chr_len, chromRadius);


       
}

function to_loci(idx, result){
    let loci_id=result[idx].loci_id;
    let loci_start=result[idx].start;
    let loci_end = result[idx].end; // should be changed
    let loci_strand = result[idx].strand;
    let gcf = result[idx].gcf;
    let nc_id = result[idx].nc_id;
    let classification = result[idx].classification;   
    let crispr_cas = {};
    if (result[idx].classification == "Medium Confidence"){
        crispr_cas = parse_5kbp(result[idx].wout_5kbp);
    }
    else if(result[idx].classification == "High Confidence"){
        kk = result[idx].win_5kbp;
        if (result[idx].wout_5kbp!="---"){
            kk = kk + "," + result[idx].wout_5kbp;
        }
        crispr_cas = parse_5kbp(kk);
    }
    var new_dict = {
        "loci_id": loci_id,
        "loci_start": loci_start,
        "loci_end" : loci_end,
        "loci_strand" : loci_strand,
        "gcf": gcf,
        "nc_id": nc_id,
        "classification": classification,
        "crispr_cas":crispr_cas,
        "details":[]
    }
    
    return new_dict;

}

function add_arrow_setting(x1, y1, strand, arrow_length, color){
    new_theta = 90;
    if(strand == "-"){
        new_theta += 180;
        
    }
    w = 6 //The width for lower part of arrow
    wa = 16//Width of arrow
    h = arrow_length // total height
    ha = 6// arrow of height
    jiancha = (wa-w)/2//Lower and higher width difference
    hs = h - ha
    var loci = s.polyline([x1, y1, x1+w, y1, x1+w, y1-hs, x1+w+jiancha, y1-hs, x+(w/2), y-h, x-jiancha, y-hs, x, y-hs]);
    loci.attr({
        fill:color,
        stroke:color,
        })
    var m = new Snap.Matrix();
    m.rotate(new_theta, x1,y1);
    loci.transform(m)
}

function add_arrow_line(s, start_x, start_y, offset, strand, color){
    x = start_x;
    y = start_y;
    l = offset;

    if(strand == "+"){
        var aca_acr = s.polyline([
            x, y-5,
            x, y+5,
            x+l-5, y+5,
            x+l-5, y+10,
            x+l,y,
            x+l-5,y-10,
            x+l-5,y-5
        ]);
        aca_acr.attr({
            fill:color,
            stroke:color,
        });
    }//if
    else{
        var aca_acr = s.polyline([
            x,y,
            x+5, y-10,
            x+5, y-5,
            x+l, y-5,
            x+l, y+5,
            x+5, y+5,
            x+5, y+10
        ]);
        aca_acr.attr({
            fill:color,
            stroke:color,
        });
    }//else
}
function add_crispr_cas(s, tag, result_json, crispr_cas_start_y, centerX, centerY, chr_len, chromRadius){
    var crispr_list={};
    var spacer_no = 1;
    var lineStart = 10000000;
    var lineEnd = -100000;
    var cas_proteins = tag[0].cas; //get cas proteins
    console.log("5kbp");
    console.log(tag);
    console.log("result json");
   console.log(result_json);
   //get crispr list, get line start and line end value.
    for(tagi in tag){
        // console.log(tag[tagi]);
        var crispr_name = tag[tagi].crispr.Name; // NC_xxx.11_1
        
        //judge whether crisprList has name
        if((crispr_name in crispr_list) == false){
            crispr_list[crispr_name] ={
                "spacer":[tag[tagi].Spacer_Pos],
                "spacer_no" : [spacer_no],
                "crispr":tag[tagi].crispr
            };
            lineStart = Math.min(tag[tagi].crispr.Start, lineStart);
            lineEnd = Math.max(tag[tagi].crispr.End, lineEnd);
            spacer_no += 1;
        }
        else{
            var s_pos = tag[tagi].Spacer_Pos;
            var temp_s_pos = crispr_list[crispr_name].spacer;
            if(temp_s_pos.includes(s_pos)){
                continue;
            }
            else{
                crispr_list[crispr_name].spacer.push(s_pos);
                crispr_list[crispr_name].spacer_no.push(spacer_no);
                spacer_no += 1;
            }
        }
    }
    var crispr_cas_pairs = [];
    var tt = [];
    console.log(crispr_list);
    let ff = [];
    for(let ci in cas_proteins){
        ff.push(cas_proteins[ci]);
    }
    for(let crispr_no in crispr_list){
        let region_start = parseInt(crispr_list[crispr_no].crispr.Start);
        let region_end = parseInt(crispr_list[crispr_no].crispr.End);
        let pi = match_crispr_cas(region_start, region_end, cas_proteins);
        if(pi!=-1){
            // cas_proteins[pi].crispr = crispr_list[crispr_no].crispr;
            crispr_cas_pairs.push({"cas":cas_proteins[pi], "crispr":crispr_list[crispr_no].crispr});
            if(!tt.includes(cas_proteins[pi])){
                tt.push(cas_proteins[pi])
                console.log(cas_proteins[pi]);
            }
        }
        
       
    }
    // console.log("hahaha");
    // console.log(crispr_cas_pairs);
    if(is_assemble==="1"){
        add_cas(s, ff, centerX, centerY, chr_len, chromRadius);
    }
    
    var fig_y = crispr_cas_start_y +50;
    draw_cas_and_crispr_array(s, crispr_cas_pairs, fig_y);
}

function draw_cas_and_crispr_array(s, crispr_cas_pairs, start_y){
    console.log("crispr_cas_pairs");
    console.log(crispr_cas_pairs);
    // console.log(start_y);
    //  cas     gap        crispr
    let y = start_y
    if(crispr_cas_pairs.length === 0){
        let error_text = "The distance between CRISPR array and the neighboring Cas locus is larger than 10kb so we won't show it!"
        let ttx = s.text(20, y, error_text);
        s.text({
            fill:'black',
            stroke:'black'
        });
    }
    else{
        for(let pi in crispr_cas_pairs){
            let cas = crispr_cas_pairs[pi].cas;
            let crispr = crispr_cas_pairs[pi].crispr;
            let cas_length = parseInt(cas.end) - parseInt(cas.start);
            let crispr_length = crispr.End - crispr.Start;
            let vs = 20;
            let ve = 1300;
            let gap = 30;
            let ratio = (ve - vs - gap) / (cas_length + crispr_length);
            let cas_length_v = cas_length * ratio;
            let crispr_length_v = crispr_length * ratio;
            // for crispr
            let spacer_fig_length = crispr_length_v / crispr.Regions.length; //e.g. 6, pixel
            console.log(`${parseInt(cas.start)} ${crispr.Start} ${parseInt(cas.end)} ${crispr.End}`);
            console.log(parseInt(cas.start)<=crispr.Start && parseInt(cas.end)>=crispr.End);
            if (parseInt(cas.end) <= crispr.Start){
                draw_cas(s, vs, cas_length_v, cas, y, ratio, parseInt(cas.start));
                
                let vss = vs + cas_length_v + gap;
                let tt = s.text(vs + cas_length_v, y + 10, "......" );
                tt.attr({
                    fill:"black",
                    stroke:"black"
                });
                let tt2 = s.text(vs + cas_length_v, y + 35, "......" );
                tt2.attr({
                    fill:"black",
                    stroke:"black"
                });
                let is_draw_ruler = true;
                draw_crispr(s, vss, crispr_length_v, crispr, y, ratio, crispr.Start, is_draw_ruler);
                
            }
            else if(crispr.End <= parseInt(cas.start)){
                let is_draw_ruler = true;
                draw_crispr(s, vs, crispr_length_v, crispr, y, ratio, crispr.Start, is_draw_ruler);
                let vss = vs + crispr_length_v + gap;
                let tt = s.text(vs + crispr_length_v, y + 10, "......" );
                tt.attr({
                    fill:"black",
                    stroke:"black"
                });
                let tt2 = s.text(vs + crispr_length_v, y + 35, "......" );
                tt2.attr({
                    fill:"black",
                    stroke:"black"
                });
                draw_cas(s, vss, cas_length_v, cas, y, ratio, parseInt(cas.start));
                
            }else if(parseInt(cas.start)<=crispr.Start && parseInt(cas.end)>=crispr.End){
                let is_draw_ruler = false;
                console.log("is_draw_ruler");
                console.log(is_draw_ruler);
                draw_cas(s, vs, cas_length_v, cas, y, ratio, parseInt(cas.start));
                let crispr_vs  = (crispr.Start - parseInt(cas.start)) * ratio + vs;
                draw_crispr(s, crispr_vs, crispr_length_v, crispr, y, ratio, crispr.Start, is_draw_ruler);
                
            }
                y = y + 100;
        }
    }//else

}
function draw_crispr(s, vs, len, crispr, y, ratio, real_start, is_draw_ruler){
    let colors = [
        "#e8d1d1",
        "#e6c1c1",
        "#b3ecf2",
        "#c5d8ed",
        "#d1beeb",
        "#fc81bf"
    ]
    let one_len = len / crispr.Regions.length;
    let next_start = vs;
    let crispr_array = crispr.Regions;
    let j = 0;
    let spacer_color = colors[j%colors.length];
    for(let i in crispr_array){
            if(["DR", "Spacer"].includes(crispr_array[i].Type)){
                var name = crispr_array[i].Type;
                if(name==="Spacer"){
                    spacer_color = colors[j%colors.length];
                    j++;
                }
                next_start = add_crispr(s, next_start, y, one_len, name, spacer_color);
            }
        }
    let cas_text = s.text(vs, y-20,`CRISPR (${crispr.Name}:${crispr.Start}-${crispr.End})`);
    cas_text.attr({
        fill:"black",
        stroke:"black",
        fontSize:10,
        strokeWidth: 0.3
        // fontWeight:"normal"
       
    });
    if (is_draw_ruler===true){
        draw_ruler(s, y+10, vs, vs+len, ratio, real_start);
    }
    
}
function draw_cas(s, vs, len, cas, y, ratio, real_start){
        // console.log(cas);
        let cas_name = cas.id;
        let start = parseInt(cas.start);
        let end = parseInt(cas.end);
        let strand = cas.strand;
        //draw arrow
        let genes = cas.Genes;
        let p_off = [10, 15, 20, 15]
        let n_off = [-10, -15, -20, -15]
        let pii = 0;
        let nii = 0;
        let mmid = (vs + len) / 2;
        var big_cas_text = s.text(mmid, y-20, `${cas_name}[${start}, ${end}]`);
        big_cas_text.attr({
            stroke:"black",
            fontSize:12,
            strokeWidth: 0.5
        });
        for( let gi in genes ){
            let arr = genes[gi];
            let name = arr.Sub_type.split("_")[0];
            // console.log(name);
            let rstart = arr.Start;
            let rend = arr.End;
            let strand = arr.Orientation;
            let vstart = (rstart - real_start) * ratio + vs;
            let vend = (rend - real_start) * ratio + vs;
            let vlen = vend - vstart;
            add_arrow_line(s, vstart, y, vlen, strand, "pink"); // vs-> v start   y const, length, strand, color
            //add text
            let off = -10;
            let off_line = 0;
            let mid = (vstart + vend)/2;
            var cas_text = s.text(mid, y+off,`${name}`);
            cas_text.attr({
                stroke:"black",
                fontSize:12,
                strokeWidth: 0.5
            });
            //draw line
            var cas_line = s.line(mid, y+off_line, mid, y+off);
            cas_line.attr({
                stroke:"black",
                fill:"black"
            });
        }
        
        //draw text
        
        //draw arrow
        draw_ruler(s, y+5, vs, vs+len, ratio, real_start)

}
function match_crispr_cas(region_start, region_end, cas_proteins){
    var closest_pi = -1;
    let least_distance = 10000;
    console.log(cas_proteins);
    console.log(`${region_start} ${region_end}`);
    for(let pi in cas_proteins){
        let cstart = parseInt(cas_proteins[pi].start);
        let cend = parseInt(cas_proteins[pi].end);
        //cas is on  the left
        if(region_start>=cend){
            let temp_dis = region_start - cend;
            if(temp_dis < least_distance){
                closest_pi = pi;
                least_distance = temp_dis;
            }
        }
        //cas if one the right
        else if(cstart >= region_end){
            let temp_dis = cstart - region_end;
            if(temp_dis < least_distance){
                closest_pi = pi;
                least_distance = temp_dis;
            }
        }
        // crispr is in the cas
        else if(cstart<=region_start && cend>=region_end){
            let temp_dis = 0;
            if(temp_dis < least_distance){
                closest_pi = pi;
                least_distance = temp_dis;
            }
        }
        else{
            console.log("error");
        }
    }
    
    return closest_pi;
}
function draw_ruler(s, line_y, line_start, line_end, ratio, real_start){
    // draw ruler
    var ruler_y = line_y + 25;
    var line  = s.line(line_start, ruler_y, line_end, ruler_y);
    // the line length is 400, acutal length is loci_length
    line.attr({
        fill: "#fff", //color
        stroke: "#000", // color
        strokeWidth: 2//width of line
    });
    var interval = 100;
    var real_interval = interval / ratio;
    // var real_start = parseInt(loci_start);
    var line_start_copy = line_start;
    while(line_start_copy <line_end){
        //draw line at line end
        var line1  = s.line(line_start_copy, ruler_y, line_start_copy, ruler_y+5); // draw kedu
        line1.attr({
            fill: "#fff", //color
            stroke: "#000", // color
            strokeWidth: 2//width of line
        });
        
        var line_text = s.text(line_start_copy, ruler_y + 16, Math.floor(real_start).toString());
        line_text.attr({
            fill:"black",
            fontSize:10
        });
        line_start_copy += interval; 
        real_start += real_interval;
    }//while
    var line1  = s.line(line_end, ruler_y, line_end, ruler_y+5); // draw kedu
    line1.attr({
        fill: "#fff", //color
        stroke: "#000", // color
        strokeWidth: 2//width of line
    });
}
          
function add_crispr(s, fig_crispr_s, fig_y, fig_width, name, spacer_color){
    
    if(name=="Spacer"){
        let c = s.rect(fig_crispr_s, fig_y-fig_width, fig_width, 2*fig_width);
        let fill = spacer_color;
        c.attr({
            "fill":fill,
            "stroke":"black"
        });
    }
    else if(name =="DR"){
        let mid = fig_width/2;
        let x1 = fig_crispr_s;
        let y1 = fig_y;
        var diamond = s.polyline(
            x1, y1,
            x1+mid, y1-fig_width,
            x1+2*mid, y1,
            x1+mid, y1+fig_width
        ); 
        diamond.attr({
            "fill":"black",
            "stroke":"black"
        });
    }
    return fig_crispr_s+fig_width;
}

function draw_aca_acr_details(s, arr2, chromRadius, centerY){
    // draw aca_acr loci lines
    var largestOffset = 20;
    for(var jj=0; jj<arr2.length; jj++){
        classes  =  arr2[jj].classification;
        if((classes == "High Confidence" || classes == "Medium Confidence") && is_assemble==="1"){
            largestOffset = 300;
        }
    }
    var offset = chromRadius + largestOffset;
    var ratio = 0.5;
    var color_aca = "#6D8764";
    var color_acr = "#716775";    
    for(i=0;i<arr2.length; i++){
        //get loci info
        loci_start = arr2[i].loci_start;//"237406"
        loci_end = arr2[i].loci_end;// "238109"
        let refseq_id = arr2[i].nc_id;
        
        if(arr2[i].details.length>7){
            ratio = 0.2;
            fontsize = 10;
        }else{
            ratio = 0.5;
            fontsize = 15;
        }
        loci_length = (arr2[i].loci_end - arr2[i].loci_start)*ratio;
        strand = arr2[i].loci_strand;
        one_crispr_cas = arr2[i].crispr_cas;
        offset = offset + 70;
        loci_id = arr2[i].loci_id;
        line_y = centerY+offset;
        my_text = `${refseq_id} Acr-Aca ${loci_id}`;
        text_len = my_text.length * 10; 
        var tmp_text = s.text(20, line_y, `${refseq_id} Acr-Aca ${loci_id}`);
        tmp_text.attr({
            fill:"black"
        });
        line_start =20 + text_len;
        line_end = line_start + loci_length; 
        var line  = s.line(line_start, line_y, line_end, line_y);
        // the line length is 400, acutal length is loci_length
        line.attr({
        fill: "#fff", //color
        stroke: "#000", // color
        strokeWidth: 2//width of line
        });
        //draw text for line
        
        //draw_aca_acr details
        text_s = line_end;
        tag = 1
        offset_obj = {"0":20, "1":-10};
        text_offset = offset_obj[tag];
        var last_a_name = "12";
        var last_arrow_length = 20;
        
        for(j=0; j<arr2[i].details.length; j++){
            var details = arr2[i].details[j];
            acr_aca = details.acr_aca;
            protein_id = details.protein_id;
            pro_start = details.start;
            pro_end = details.end;
            fig_start = (pro_start - loci_start)*ratio + line_start;
            fig_end = (pro_end - loci_start)*ratio + line_start;
            if (acr_aca == "Acr"){
                c = color_acr
                a_name = acr_aca + ":" + details.protein_id;
            }
            else{
                c = color_aca;
                a_name = "Aca" + ":" + details.protein_id;
            }
            //Every character length is 6.47. If last aca/acr exceed the length of word, we use xor for the current ID. 
            if(j>0){
                var word_length = last_a_name.length * 6.47;
                if(word_length> last_arrow_length){
                    tag = tag ^ 1;
                    text_offset = offset_obj[tag];
                }
            }
            arrow_length = parseInt(fig_end - fig_start);
            add_arrow_line(s, fig_start, line_y, arrow_length, strand, c);
            var tmp_text = s.text(fig_start, line_y + text_offset, a_name);
            tmp_text.attr({
                fill:c,
                fontSize:12
            });
            //let next loop to use the following variables in order to compare
            last_a_name = a_name;
            last_arrow_length = arrow_length;
        }
        // draw ruler
        var ruler_y = line_y + 25;
        var line  = s.line(line_start, ruler_y, line_end, ruler_y);
        // the line length is 400, acutal length is loci_length
        line.attr({
            fill: "#fff", //color
            stroke: "#000", // color
            strokeWidth: 2//width of line
        });
        var interval = 100;
        var real_interval = interval / ratio;
        var real_start = parseInt(loci_start) - contigs[refseq_id].start;
        var line_start_copy = line_start;
        while(line_start_copy <line_end){
            //draw line at line end
            var line1  = s.line(line_start_copy, ruler_y, line_start_copy, ruler_y+5); // draw kedu
            line1.attr({
                fill: "#fff", //color
                stroke: "#000", // color
                strokeWidth: 2//width of line
            });
            
            var line_text = s.text(line_start_copy-20, ruler_y + 16, Math.floor(real_start).toString());
            line_text.attr({
                fill:"black",
                fontSize:10
            });
            line_start_copy += interval; 
            real_start += real_interval;
        }//while
        var line1  = s.line(line_end, ruler_y, line_end, ruler_y+5); // draw kedu
        line1.attr({
            fill: "#fff", //color
            stroke: "#000", // color
            strokeWidth: 2//width of line
            });
    }//for acr_aca_loci


}



function get_arr_from_result(result){
    let arr=[];
	let temp_loci = to_loci(0, result);	
	for(let qq=0; qq<result.length; qq++){
        if (temp_loci.loci_id != result[qq].loci_id){
            arr.push(temp_loci);
			temp_loci = to_loci(qq, result);
		}
        temp_loci.loci_end = result[qq].end;
		temp_loci.details.push(result[qq]);
    }
    arr.push(temp_loci);
    return arr;
}


function draw_arc(s, contigs, centerX, centerY, chr_len, chromRadius){
    var arc_color = ["red", "blue", "orange"];
    let i = 0;
    for (let refseq in contigs){
        let start_pos = contigs[refseq].start;
        let end_pos = contigs[refseq].end;
        [sxx, syy, theta1] = get_locus_xy(start_pos, chr_len, centerX, centerY, chromRadius);
        [exx, eyy, theta2] = get_locus_xy(end_pos, chr_len, centerX, centerY, chromRadius);
        sweep_flag="1";
        bigArcFlag = "1";
        let mypath = `M${sxx} ${syy} A${chromRadius} ${chromRadius}, 0, ${bigArcFlag} ${sweep_flag}, ${exx} ${eyy}`;
        console.log(mypath);
        var arc_fig  =  s.path(mypath);
        arc_fig.attr({
            stroke:"blue",
            fill:"none"
        });
        i++;
    }
}
function add_location_to_loci(arr, chr_len, centerX, centerY, chromRadius){
	// draw Aca-Acr Loci using arr
    
    //word length is 200.57
    for(let i=0;i<arr.length;i++){
        if(i==0){
            [xx,yy, new_theta]= get_locus_xy(arr[i].loci_start, chr_len, centerX, centerY, chromRadius);
            arr[i]['fig']={
                'x1':xx,
                'y1':yy,
                'new_theta':new_theta,
                'text_x1':xx,
                'text_y1':yy
            }
        }
        else{
            [xx,yy, new_theta]= get_locus_xy(arr[i].loci_start, chr_len, centerX, centerY, chromRadius);
            arr[i]['fig']={
                'x1':xx,
                'y1':yy,
            }
            var text_y2 = arr[i-1]['fig'].text_y1;//ori
            var theta2 = arr[i-1]['fig'].new_theta;//ori
            var text_x1 = xx;
            var text_y1 = yy;
            var theta1 = new_theta;
            //adjust text y
            if((theta2>180 && theta1>180 || theta2<180 && theta1<180)){
                text_y1 = get_cur_texty(text_y2, text_y1, theta2, theta1,20);
            }
            //adjust text x
            //Circile corrdinate: (x - x0)^2 + (y - y0)^2 = r^2
            var delta_y = text_y1 - centerY;
            var delta_x = text_x1 - centerX;
            var delta_y_square = Math.pow(delta_y, 2);
            if(delta_y_square < Math.pow(chromRadius, 2)){
                var text_radius = Math.sqrt(Math.pow(chromRadius, 2) - delta_y_square); 
                if(delta_x > 0){
                    text_x1 = centerX + text_radius;
                }
                else{
                    text_x1 = centerX - text_radius;
                    text_x1 -= 200;
                }            
            }
            else{ 
                if(delta_x < 0){
                    text_x1 -= 200;
                } 
            }
            
            arr[i]['fig']['new_theta'] = new_theta;
            arr[i]['fig']['text_x1'] = text_x1;
            arr[i]['fig']['text_y1'] = text_y1;
        }//else
    
    }
    return arr;
}
function draw_graphic(result2, s, chromRadius, is_assemble, centerX, centerY, crispr_cas_start_y){
    //1. draw circular chromosome
    if(is_assemble==="1"){
        let posCircle = s.circle(centerX, centerY, chromRadius); // draw chrom using circle
        posCircle.attr({
        fill: "#fff", //Circle color
        stroke: "#000", // Circle line color
        strokeWidth: 2 //width of line
        });
         //draw original point
        let text_ori_point = s.text(centerX, centerY+chromRadius+30, `Start/End point:0/${chr_len}`);
        text_ori_point.attr({
            fill: "#000", //color
            stroke: "#000", // color
            strokeWidth: 0.5
        })
    }
    //2. draw arc 
    
    let arr2 = get_arr_from_result(result2);
    if(is_assemble==="1"){
        arr2 = add_location_to_loci(arr2, chr_len, centerX, centerY, chromRadius);
    }
	var drawn=false;//crispr_cas has drawn or not
    draw_aca_acr_details(s, arr2, chromRadius, centerY);
    // var crispr_cas_tag;
	for(let arri=0; arri<arr2.length; arri++){
        tag = arr2[arri];
        if(is_assemble==="1"){
            draw_loci(s, tag);
        }
		if(!jQuery.isEmptyObject(tag.crispr_cas) && drawn === false && json_content!="none"){
            add_strand_crispr(s, tag.crispr_cas, json_content, centerX, centerY, chr_len, chromRadius, crispr_cas_start_y);
            drawn = true;
		}
    }
    

}

