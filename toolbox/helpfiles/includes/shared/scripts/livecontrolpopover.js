$(function () {
    $('[data-toggle="popover"]').popover({
        title: 'Interactive Control',
        content: function() {
            var cmd = $(this).attr('data-examplename');
            var content = '<span>In live scripts, controls like this one let you set the value of a variable interactively. To use the controls in this example, <a class="no-matlab" href="matlab:' + cmd + '">open the live script</a>.</span>';
            return content;
        }
    })
    $(window).trigger('popover_added');
});
