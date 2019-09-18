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
        .attr("class", "form-group col-3")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control required-form-control file-ctrl-group py-1")
        .attr("id", "file-ctrl-group")
        .attr("name", "file-ctrl-group")
        .attr("placeholder", "Name of control group")
        .attr("required", true);
    newFormRow.append('div')
        .attr("class", "form-group col-3")
        .append('input')
        .attr("type", "text")
        .attr("class", "form-control required-form-control file-ref-gene py-1")
        .attr("id", "file-ref-gene")
        .attr("name", "file-ref-gene")
        .attr("placeholder", "Name of reference gene")
        .attr("required", true);
    let newSel = newFormRow.append('div')
        .attr("class", "form-group col-1")
        .append('select')
        .attr("class", "form-control required-form-control file-repeat py-1")
        .attr("id", "file-repeat")
        .attr("name", "file-repeat")
        .attr("required", true);
    newSel.append('option')
        .attr("value", "")
        .attr("selected", true)
        .attr("hidden", true)
        .text("2 or 3");        
    newSel.append('option')
        .attr("value", "2")
        .text("2");         
    newSel.append('option')
        .attr("value", "3")
        .text("3");
    newFormRow.append('div')
        .attr("class", "form-group col-1")
        .append('i')
        .attr("class", "fas fa-minus-circle mt-2 pt-1");

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
    if (d3.select(".params-bar-colors").node() != null) {
        d3.select(".params-bar-colors").remove();
    }
    
    // Insert 'div' element above "Bar sorting" input box
    let sampleGroupColors = d3.select(".upload")
        .insert('div', ".params-t-sort")
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

// Event listener when changing "Specific sorting on sample names"
d3.select(".params-s-sort").selectAll('input').on("change", function() {    
    
    // Variable for the selected checkbox
    let checkboxSel = d3.select(this).node().value;
    
    // Check if "specific sorting on sample names" is selected
    if (checkboxSel == "yes") {

        // Make sure the opposite checkbox is unchecked
        d3.select("#params-s-sort-no")
            .property("checked", false);        

        // Specific sorting on sample names is not pre-checked
        if (d3.select(".params-s-sort-range").node() === null) {
            // Insert 'div' element above "Total bar thold for h-plot" input box
            let sampleSortRange = d3.select(".upload")
                .insert('div', ".params-total-bars")
                .attr("class", "form-row plot-params params-s-sort-range px-2")
                .append('div')
                .attr("class", "form-group col-10 d-flex");                
            // Append position range on sample name for sorting
            sampleSortRange.append('label')
                .attr("for", "s-sort-start")
                .attr("class", "px-1 mb-3")
                .text("Range (couting from 0): ");          
            sampleSortRange.append('span')
                .html("&emsp;");            
            sampleSortRange.append('input')
                .attr("type", "number")
                .attr("id", "s-sort-start")
                .attr("class", "form-control col-2")
                .attr("name", "s-sort-start")
                .attr("min", "-10000000")
                .attr("step", "1")
                .attr("placeholder", "Starting position on sample name for sorting");
            sampleSortRange.append('span')
                .html("&ensp;");         
            sampleSortRange.append('input')
                .attr("type", "number")
                .attr("id", "s-sort-end")
                .attr("class", "form-control col-2")
                .attr("name", "s-sort-end")
                .attr("min", "-10000000")
                .attr("step", "1")
                .attr("placeholder", "Ending position on sample name for sorting");                            
        // Specific sorting on sample names is pre-checked
        } else {
            // Remove element for sample sorting range
            d3.select(".params-s-sort-range").remove();
        }
    
    // Not "specific sorting on sample names" is selected
    } else {

        // Make sure the opposite checkbox is unchecked
        d3.select("#params-s-sort-yes")
            .property("checked", false);
        // Remove element for sample sorting range if exists
        if (d3.select(".params-s-sort-range").node() != null) {
            d3.select(".params-s-sort-range").remove();
        }        
         
    }

});

// Make 'img' visible after clicking "PLOT" button
d3.select("#form-submit-plot").on("click", function() {
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