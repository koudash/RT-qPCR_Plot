// Event listener when clicking "Add data" button
d3.select("#file-add").on("click", function() {

    // Prevent self-refreshing
    d3.event.preventDefault();

    // Insert new data row above "Add data" button
    let newFormRow = d3.select(".upload")
        .insert('div', ".but-bottom")
        .attr("class", "form-row upload-row px-2");
    newFormRow.append('div')
        .attr("class", "form-group col-4")
        .append('input')
        .attr("type", "file")
        .attr("class", "form-control required-form-control file-sel py-1")
        .attr("id", "file-sel")
        .attr("name", "file-sel")
        .attr("accept", ".csv")
        .attr("required", true);
    newFormRow.append('div')
        .attr("class", "form-group col-4")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control required-form-control file-ref-gene py-1")
        .attr("id", "file-ref-gene")
        .attr("name", "file-ref-gene")
        .attr("placeholder", "Name of reference gene")
        .attr("required", true);
    newFormRow.append('div')
        .attr("class", "form-group col-4")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control required-form-control file-ctrl-group py-1")
        .attr("id", "file-ctrl-group")
        .attr("name", "file-ctrl-group")
        .attr("placeholder", "Name of control group")
        .attr("required", true);
    
        // Note that event listeners for row deletion and text color change have to be updated once new row is added
        // Delete row once "delete" icon is clicked
        d3.selectAll(".fa-minus-circle").on("click", function() {
            $(this).parents()[1].remove();
        });
        // Change color of texts in required input box when filled in
        d3.selectAll(".required-form-control").on("change", function() {
            this.style.color = "#5CB85C";
        });
        
});

// Event listener when changing "No. of sample groups"
d3.select("#params-sample-groups").on("change", function() {

    // Remove pre-existing selections for bar colors should there be any
    if (d3.select(".params-bar-colors")) {
        d3.select(".params-bar-colors").remove();
    }
    
    // Insert 'div' element above "Bar sorting" input box
    let sampleGroupColors = d3.select(".upload")
        .insert('div', ".params-sort")
        .attr("class", "form-row plot-params params-bar-colors px-2");

    // Variable for the text in "No. of sample groups" box
    // https://stackoverflow.com/questions/17527872/get-value-of-input-element-in-event-listener-with-d3
    let sampleGroupsNum = d3.select(this).node().value;
    
    // Append selections for bar colors for each sample group
    for (i=0; i<sampleGroupsNum; i++) {
        sampleGroupColors.append('label')
            .attr("for", `params-bar-color-${i}`)
            .attr("class", "individual-bar-color px-1 mb-3")
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

    // Sorting by value is available if there are only two sample groups (one ctrl and one treatment groups)
    if (sampleGroupsNum == 2) {
        d3.selectAll(".sort-values").attr("class", "sort-values shown");
    } else {
        d3.selectAll(".sort-values").attr("class", "sort-values hidden");
    }

});

// Make 'img' visible after clicking "PLOT" button
d3.select("#form-submit").on("click", function() {
    d3.select('.bar-graph')
        .attr("class", "bar-graph shown");
});

// Event listener for optional inputs
$(".form-subheader-2").click(function() {
    // Toggling to decompress/collapse optional inputs
    // https://stackoverflow.com/questions/13098863/fold-div-on-click
    $(".optional-form-group").toggle(500);
    // Toggling between "up" and "down" icons
    // https://stackoverflow.com/questions/3965479/trying-to-toggle-between-two-icons
    $(".fa-angle-double-down").toggleClass("fa-angle-double-up");
});

// Change color of texts in required input box when filled in
d3.selectAll(".required-form-control").on("change", function() {
    this.style.color = "#5CB85C";
});