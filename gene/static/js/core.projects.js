$(document).ready(function() {
    var changeProject = "";
    var changeSubsamplingAnchor = null;
    hideLoading();

    // Popovers on the create page
    $('[data-toggle="popover"]').popover();

    $(".download").click(function() {
        var type = $(this).data("downloadtype");
        var project = $(this).data("project");
        console.log("Download file " + type + " in project " + project);
        window.open("download?pid=" + project + "&type=" + type, '_blank');
    });

    $(".add-example-project").click(function() {
        $.ajax({
            type: "POST",
            url: "/loadExampleProject",
            success: function(result) {
                var baseURL = '//' + location.host + location.pathname;
                window.location.href = baseURL;
            },
            error: function(err) {
                console.log(err);
            }
        });
    });

    $(".project-delete").click(function() {
        var project_name = $(this).data("projectname");
        var project = $(this).data("project");
        bootbox.confirm(
            "Are you sure you want to delete " + project_name + "?",
            function(result) {
                if (result) {
                    var data = {
                        project: project,
                        delete: "delete"
                    };

                    $.ajax({
                        type: "POST",
                        url: "/deleteProject",
                        data: data,
                        success: function(result) {
                            // TODO: Fix to use IDs
                            $("#p-" + project).remove();
                        },
                        error: function(err) {
                            console.log(err);
                        }
                    });
                }
            }
        );
    });

    $(".project-trig-otu").click(function() {
        var project = $(this).data("project");
        $(this)
            .siblings(".project-replace-otu")
            .trigger("click");
    });
    $(".project-trig-metadata").click(function() {
        var project = $(this).data("project");
        $(this)
            .siblings(".project-replace-metadata")
            .trigger("click");
    });

    $(".project-replace-otu").change(function() {
        upload($(this));
    });
    $(".project-replace-metadata").change(function() {
        upload($(this));
    });

    $("#change-close").click(function() {
        $("#change-box").hide();
        $("#blackout").hide();
    });

    $(".change-cancel").click(function() {
        $("#change-box").hide();
        $("#blackout").hide();
    });
});
