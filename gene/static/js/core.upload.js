$(document).ready(function() {
    // Hides loading bar
    $("#loading").hide();

    // Global variables
    var isNameValid = false;
    var otuTableUploaded = false;
    var otuTaxonomyMappingUploaded = false;
    var otuMetadataUploaded = false;

    // Event listener for main create button
    $("#upload-submit").click(function() {
        $("#loading").show();
    });

    // The event listener for the file upload
    document
        .getElementById("inputName")
        .addEventListener("change", nameChange, false);
    document
        .getElementById("otuTable")
        .addEventListener("change", uploadOTUTable, false);
    document
        .getElementById("otuMetadata")
        .addEventListener("change", uploadMetadata, false);

    function nameChange() {
        $("#projectName").val($("#inputName").val());
        checkComplete();
    }

    function checkComplete() {
        if ($("#projectName").val() !== "") {
            isNameValid = true;
        } else {
            isNameValid = false;
        }

        if (
            isNameValid &&
            (otuTableUploaded && otuMetadataUploaded)
        ) {
            $("#upload-submit").prop("disabled", false);
        } else {
            $("#upload-submit").prop("disabled", true);
        }
    }

    function uploadOTUTable() {
        upload("otuTableForm");
    }

    function uploadMetadata() {
        upload("otuMetadataForm");
    }

    // Catch the form submit and upload the files
    function upload(formID) {
        var form = $("#" + formID)[0];
        var filename = "";
        var formData = new FormData(form);

        if (formID == "otuTableForm") {
            $("#otuTableLoading").show();
            filename = $("#otuTable").val();
        }
        if (formID == "otuMetadataForm") {
            $("#otuMetadataLoading").show();
            filename = $("#otuMetadata").val();
        }

        filename = filename.split(/(\\|\/)/g).pop();

        $.ajax({
            url: "/upload",
            type: "POST",
            data: formData,
            cache: false,
            processData: false, // Don't process the files
            contentType: false, // Set content type to false as jQuery will tell the server its a query string request
            success: function(data, textStatus, jqXHR) {
                console.log("Success");
                if (formID == "otuTableForm") {
                    $("#otuTableLoading").hide();
                    $("#otuTableOK").show();
                    $("#otuTableText").text("Replace");
                    $("#projectOTUTableName").val(filename);
                    otuTableUploaded = true;
                }
                if (formID == "otuMetadataForm") {
                    $("#otuMetadataLoading").hide();
                    $("#otuMetadataOK").show();
                    $("#otuMetadataText").text("Replace");
                    $("#projectSampleIDName").val(filename);
                    otuMetadataUploaded = true;
                }

                checkComplete();
            },
            error: function(jqXHR, textStatus, errorThrown) {
                console.log("Error: " + textStatus);
                alert("An error occurred");
            }
        });
    }
});
