
// // Variable for 'button' element to add new file record
// var butAddFile = d3.select("#file-add");

d3.select("#file-add").on("click", function() {

    // Prevent self-refreshing
    d3.event.preventDefault();

    let newFormRow = d3.select(".upload")
        .insert('div', ".but-bottom")
        .attr("class", "form-row upload-row");

    newFormRow.append('div')
        .attr("class", "form-group col-5")
        .append('input')
        .attr("type", "file")
        .attr("class", "form-control file-sel py-1")
        .attr("name", "file-sel")
        .attr("accept", ".csv");
        
    newFormRow.append('div')
        .attr("class", "form-group col-3")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control fiel-ref-gene py-1")
        .attr("name", "file-ref-gene")
        .attr("placeholder", "Name of reference gene");

    newFormRow.append('div')
        .attr("class", "form-group col-3")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control fiel-ctrl-group py-1")
        .attr("name", "file-ctrl-group")
        .attr("placeholder", "Name of control group");        

    newFormRow.append('div')
        .attr("class", "form-group col-1")
        .append('i')
        .attr("class", "fas fa-minus-circle mt-2 pt-1");
    
        // Note that row deletion event listener has to be updated once new row is added
        d3.selectAll(".fa-minus-circle").on("click", function() {

            $(this).parents()[1].remove();

        });
        
});

d3.select("#params-sample-groups").on("change", function() {

    if (d3.select(".params-bar-colors")) {
        d3.select(".params-bar-colors").remove();
    }
    
    let sampleGroupColors = d3.select(".upload")
        .insert('div', ".params-sort")
        .attr("class", "form-row plot-params params-bar-colors");

    // https://stackoverflow.com/questions/17527872/get-value-of-input-element-in-event-listener-with-d3
    let sampleGroupsNum = d3.select(this).node().value;
    
    for (i=0; i<sampleGroupsNum; i++) {

        sampleGroupColors.append('label')
            .attr("for", `params-bar-color-${i}`)
            .attr("class", "px-1 mb-3")
            .text(`Bar color-${i}`);

        sampleGroupColors.append('span')
            .html(":&ensp;");            

        sampleGroupColors.append('input')
            .attr("type", "color")
            .attr("id", `params-bar-color-${i}`)
            .attr("name", "params-bar-colors");

        sampleGroupColors.append('span')
            .html("&emsp;");

    }

    if (sampleGroupsNum == 2) {

       

        d3.selectAll(".sort-values").attr("class", "sort-values shown");
    } else {
       
        d3.selectAll(".sort-values").attr("class", "sort-values hidden");
    }




});

// Make 'img' visible after data and plotting params have been submitted
d3.select("#form-submit").on("click", function() {

    d3.select('.bar-graph')
        .attr("class", "bar-graph shown");

});
