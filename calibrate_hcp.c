/*
 *  @(#) $Id: calibrate_hcp.c 2014-05-08 $
 *  Copyright (C) 2014 Jeffrey J. Schwartz.
 *  E-mail: schwartz@physics.ucla.edu
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
 */

/*
 *  This program serves to calibrate SPM images using the 2D FFT
 *  of a known hexagonal lattice.  Selecting two peaks in the FFT
 *  allows one to determine the scaling factor by which to stretch or
 *  squeeze the image so as to to make the measured lattice conform
 *  to the known lattice parameters.
 */

#include "config.h"
#include <gtk/gtk.h>
#include <app/gwyapp.h>
#include <app/gwymoduleutils.h>
#include <libprocess/stats.h>
#include <libprocess/filters.h>
#include <libprocess/inttrans.h>
#include <libprocess/datafield.h>
#include <libgwyddion/gwymath.h>
#include <libgwydgets/gwydataview.h>
#include <libgwydgets/gwydgetutils.h>
#include <libgwydgets/gwynullstore.h>
#include <libgwydgets/gwylayer-basic.h>
#include <libgwydgets/gwyradiobuttons.h>
#include <libgwymodule/gwymodule-process.h>

#define CALIBRATE_HCP_RUN_MODES (GWY_RUN_INTERACTIVE)

typedef struct _GwyToolLevel3      GwyToolLevel3;

typedef struct _GwyToolLevel3Class GwyToolLevel3Class;

struct _GwyToolLevel3
{
    GwyPlainTool parent_instance;
    GtkTreeView *treeview;
    GtkTreeModel *model;
    GtkObject *radius;
    gint32 rpx;
    GtkWidget *instant_apply;
    GtkWidget *set_zero;
    GtkWidget *apply;
    GType layer_type_point;
};

enum
{
    COLUMN_I,
    COLUMN_X,
    COLUMN_Y,
    COLUMN_Z,
    NCOLUMNS
};

enum
{
    PREVIEW_SIZE = 512
};

typedef enum {
    ZOOM_1 = 1,
    ZOOM_2 = 2,
} ZoomMode;

typedef struct {
    gdouble lower;
    gdouble upper;
    gdouble lattice;
    gdouble Xscale;
    gdouble Yscale;
    gboolean Xwarning;
    gboolean Ywarning;
    ZoomMode zoom_mode;
} ThresholdArgs;

typedef struct {
    gdouble min, max;
} ThresholdRanges;

typedef struct {
    ThresholdArgs *args;
    const ThresholdRanges *ranges;
    GtkWidget *dialog;
    GtkWidget *view;
    GtkWidget *lower;
    GtkWidget *upper;
    GtkWidget *xscale;
    GtkWidget *yscale;
    GtkWidget *xwarning;
    GtkWidget *ywarning;
    GtkWidget *warning;
    GwyContainer *mydata;
    GwyContainer *container;
    GwyDataField *ofield;
    GwyDataField *offt;
    GwyDataField *disp_data;
    GwyDataField *dfield;
    gint id;
    GwySelection *selection;
    GwySIValueFormat *original_XY_Format;
    GwySIValueFormat *XY_Format;
    GwySIValueFormat *Z_Format;
    GwyToolLevel3 *tool;
    GtkWidget *lattice;
    gdouble p[2][3];
    GSList *zoom_mode_radios;
} ThresholdControls;

static gboolean module_register             (void);

static void     calibrate_hcp               (GwyContainer *data, GwyRunType run);
static void     perform_fft                 (GwyDataField *dfield,
                                                GwyContainer *data);
static void     selection_changed           (ThresholdControls *controls);
static void     clear_points                (ThresholdControls *controls);
static void     peak_find                   (ThresholdControls *controls,
                                                gdouble *point, guint idx);
static void     calibrate_update_scales     (ThresholdControls *controls);
static void     calibration_get_factors     (ThresholdControls *controls);
static void     check_warnings              (ThresholdControls *controls);
static void     calibrate_do                (ThresholdControls *controls);
static void     calibrate_create_output     (GwyContainer *data, 
                                                GwyDataField *dfield,
                                                ThresholdControls *controls);
static void     calibrate_hcp_dialog        (ThresholdArgs *args,
                                                ThresholdRanges *ranges,
                                                GwyContainer *data,
                                                GwyDataField *dfield,
                                                gint id, GwyToolLevel3 *tool);
static void     threshold_set_to_full_range(ThresholdControls *controls);
static void     threshold_lower_changed    (ThresholdControls *controls);
static void     threshold_upper_changed    (ThresholdControls *controls);
static void     preview                    (ThresholdControls *controls);
static void     threshold_do               (const ThresholdArgs *args,
                                            GwyDataField *dfield);
static void     threshold_load_args        (GwyContainer *settings,
                                                ThresholdArgs *args,
                                                GwyToolLevel3 *tool);
static void     threshold_save_args        (GwyContainer *settings,
                                                ThresholdArgs *args,
                                                GwyToolLevel3 *tool);
static void     scale_entry_attach         (ThresholdControls *controls,
                                                GtkTable *table, gint row);
static void     threshold_lattice_changed  (ThresholdControls *controls);
static void     fft_postprocess            (GwyDataField *dfield);
static void     set_dfield_modulus         (GwyDataField *re, GwyDataField *im,
                                                GwyDataField *target);
static void     zoom_mode_changed          (GtkToggleButton *button,
                                                ThresholdControls *controls);
static void     xscale_changed             (ThresholdControls *controls);
static void     yscale_changed             (ThresholdControls *controls);
static void     gwy_tool_level3_render_cell(GtkCellLayout *layout,
                                                GtkCellRenderer *renderer,
                                                GtkTreeModel *model,
                                                GtkTreeIter *iter,
                                                gpointer user_data);
static void gwy_tool_level3_radius_changed (GwyToolLevel3 *tool);
static void radio_buttons_attach_to_table  (GSList *group,
                                                GtkTable *table, gint row);

static const ThresholdArgs threshold_defaults = {
    0.0, 0.0, 0.000000001, 1.0, 1.0, FALSE, FALSE, 1
};

static GwyModuleInfo module_info = {
    GWY_MODULE_ABI_VERSION, &module_register,
    N_("Tool to calibrate and adjust the lateral dimensions of a scanning "
    "probe microscope image using a known hexagonal close-packed lattice."),
    "Jeffrey J. Schwartz <schwartz@physics.ucla.edu>",
    "1.0",
    "Jeffrey J. Schwartz",
    "May 2014",
};

GWY_MODULE_QUERY(module_info)

static gboolean
module_register(void)
{
    gwy_process_func_register("calibrate_hcp",
                (GwyProcessFunc)&calibrate_hcp,
                N_("/_Correct Data/_Calibrate HCP"),
                NULL, CALIBRATE_HCP_RUN_MODES, GWY_MENU_FLAG_DATA,
                N_("Calibrate image against known HCP lattice"));
    return TRUE;
}

static void
calibrate_hcp(GwyContainer *data, GwyRunType run)
{
    ThresholdArgs args;
    ThresholdRanges ranges;
    GwyDataField *dfield;
    GQuark quark;
    gint id;
    GwyToolLevel3 tool;
    tool.rpx = 3;
    g_return_if_fail(run & CALIBRATE_HCP_RUN_MODES);
    threshold_load_args(gwy_app_settings_get(), &args, &tool);
    gwy_app_data_browser_get_current(GWY_APP_DATA_FIELD, &dfield,
                                     GWY_APP_DATA_FIELD_ID, &id,
                                     GWY_APP_DATA_FIELD_KEY, &quark, 0);
    g_return_if_fail(dfield);
    if (run == GWY_RUN_INTERACTIVE)
    {
        calibrate_hcp_dialog(&args, &ranges, data,
            gwy_data_field_duplicate(dfield), id, &tool);
        gwy_data_field_data_changed(dfield);
    }
}

static void
threshold_format_value(ThresholdControls *controls,
                       GtkEntry *entry, gdouble value)
{
    gchar *s;
    s = g_strdup_printf("%.*f",
                        controls->original_XY_Format->precision+1,
                        value/controls->original_XY_Format->magnitude);
    gtk_entry_set_text(GTK_ENTRY(entry), s);
    g_free(s);
}

static GtkWidget*
threshold_entry_attach(ThresholdControls *controls,
                       GtkTable *table, gint row,
                       gdouble value, const gchar *name)
{
    GtkWidget *label, *entry;
    label = gtk_label_new_with_mnemonic(name);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    entry = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(entry, TRUE);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 8);
    threshold_format_value(controls, GTK_ENTRY(entry), value);
    gtk_table_attach(table, entry, 1, 3, row, row+1, GTK_FILL, 0, 0, 0);
    label = gtk_label_new(controls->original_XY_Format->units);
    gtk_label_set_markup(GTK_LABEL(label),
                                    controls->original_XY_Format->units);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 3, 4, row, row+1, GTK_FILL, 0, 0, 0);
    return entry;
}

static void
scale_entry_attach(ThresholdControls *controls,
                       GtkTable *table, gint row)
{
    GtkWidget *label;
    label = gtk_label_new("\t\tX: ");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    controls->xscale = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(controls->xscale, FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(controls->xscale), 5);
    gtk_table_attach(table, controls->xscale, 1, 3,
                                row, row+1, GTK_FILL, 0, 0, 0);
    controls->xwarning = gtk_label_new(NULL);
    gtk_label_set_width_chars(GTK_LABEL(controls->xwarning), 2);
    gtk_misc_set_alignment(GTK_MISC(controls->xwarning), 0.0, 0.5);
    gtk_table_attach(table, controls->xwarning, 3, 4,
                                row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    label = gtk_label_new("\t\tY: ");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    controls->yscale = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(controls->yscale, FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(controls->yscale), 5);
    gtk_table_attach(table, controls->yscale, 1, 3,
                                row, row+1, GTK_FILL, 0, 0, 0);
    controls->ywarning = gtk_label_new(NULL);
    gtk_label_set_width_chars(GTK_LABEL(controls->ywarning), 2); 
    gtk_misc_set_alignment(GTK_MISC(controls->ywarning), 0.0, 0.5);
    gtk_table_attach(table, controls->ywarning, 3, 4,
                                row, row+1, GTK_FILL, 0, 0, 0);
    gwy_widget_set_activate_on_unfocus(controls->xscale, TRUE);
    gwy_widget_set_activate_on_unfocus(controls->yscale, TRUE);
    gtk_entry_set_text(GTK_ENTRY(controls->xscale), "1.0");
    gtk_entry_set_text(GTK_ENTRY(controls->yscale), "1.0");
    g_signal_connect_swapped(controls->xscale, "activate",
                            G_CALLBACK(xscale_changed), controls);
    g_signal_connect_swapped(controls->yscale, "activate",
                            G_CALLBACK(yscale_changed), controls);
}

static void
calibrate_hcp_dialog(ThresholdArgs *args, ThresholdRanges *ranges,
                 GwyContainer *data, GwyDataField *dfield,
                 gint id, GwyToolLevel3 *tool)
{
    GtkWidget *dialog, *hbox, *button, *label;
    GtkTable *table;
    GwyVectorLayer *vlayer;
    ThresholdControls controls;
    gint response;
    GwyPixmapLayer *layer;
    gint row;
    controls.ofield = gwy_data_field_duplicate(dfield);
    controls.container = data;
    controls.id = id;    
    controls.args = args;
    controls.ranges = ranges;
    controls.dfield = dfield;
    controls.tool = tool;
    controls.original_XY_Format = gwy_data_field_get_value_format_xy
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    controls.mydata = gwy_container_new();
    perform_fft(controls.dfield, controls.mydata);
    controls.offt = gwy_data_field_duplicate(controls.dfield);
    controls.disp_data = gwy_data_field_duplicate(controls.dfield);
    dfield = gwy_data_field_duplicate(controls.dfield);
    gwy_data_field_get_min_max(dfield, &ranges->min, &ranges->max);
    controls.XY_Format = gwy_data_field_get_value_format_xy
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    controls.Z_Format = gwy_data_field_get_value_format_z
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    dialog = gtk_dialog_new_with_buttons(_("Calibrate HCP"), NULL, 0,
                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_OK);
    controls.dialog = dialog;
    hbox = gtk_hbox_new(FALSE, 2);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox,
                       FALSE, FALSE, 4);
    table = GTK_TABLE(gtk_table_new(2, 1, FALSE));
    gtk_table_set_row_spacings(table, 2);
    gtk_table_set_col_spacings(table, 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);
    gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(table), TRUE, TRUE, 4);
    label = gtk_label_new("FFT of data");
    gtk_label_set_markup(GTK_LABEL(label),
                "<b>FFT of data</b>\nModulus, Hanning window, subtract mean");
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment(GTK_MISC(label), 0.5, 0.5);
    gtk_table_attach(table, label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0);
    gwy_app_sync_data_items(data, controls.mydata, id, 0, FALSE,
                GWY_DATA_ITEM_PALETTE, GWY_DATA_ITEM_MASK_COLOR,
                GWY_DATA_ITEM_RANGE, GWY_DATA_ITEM_REAL_SQUARE, 0);
    gwy_container_set_object_by_name(controls.mydata, "/0/data", dfield);
    controls.view = gwy_data_view_new(controls.mydata);
    layer = gwy_layer_basic_new();
    g_object_set(layer, "data-key", "/0/data",
                 "gradient-key", "/0/base/palette",
                 "range-type-key", "/0/base/range-type",
                 "min-max-key", "/0/base", NULL);
    gwy_data_view_set_data_prefix(GWY_DATA_VIEW(controls.view), "/0/data");
    gwy_data_view_set_base_layer(GWY_DATA_VIEW(controls.view), layer);
    gwy_set_data_preview_size(GWY_DATA_VIEW(controls.view), PREVIEW_SIZE);
    vlayer = g_object_new(g_type_from_name("GwyLayerPoint"),
                  "selection-key", "/0/select/point", NULL);
    gwy_data_view_set_top_layer(GWY_DATA_VIEW(controls.view), vlayer);
    controls.selection = gwy_vector_layer_ensure_selection(vlayer);
    gwy_selection_set_max_objects(controls.selection, 2);
    g_signal_connect_swapped(controls.selection, "changed",
                         G_CALLBACK(selection_changed), &controls);
    gtk_table_attach(table, controls.view, 0, 1, 1, 2, GTK_FILL, 0, 0, 0);
    label = gtk_label_new(
        "Select two peaks in the first hexagonal ring around center");
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment(GTK_MISC(label), 0.5, 0.5);
    gtk_table_attach(table, label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0);
    table = GTK_TABLE(gtk_table_new(5, 4, FALSE));
    gtk_table_set_row_spacings(table, 2);
    gtk_table_set_col_spacings(table, 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);
    gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(table), TRUE, TRUE, 4);
    row = 0;
    label = gtk_label_new("Display Zoom: ");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Zoom:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    row++;    
    controls.zoom_mode_radios
        = gwy_radio_buttons_createl(G_CALLBACK(zoom_mode_changed), &controls,
                                    controls.args->zoom_mode,
                                    _("×1"), ZOOM_1,
                                    _("×2"), ZOOM_2,
                                    NULL);
    radio_buttons_attach_to_table(controls.zoom_mode_radios, table, row);
    row++;
    label = gtk_label_new("Specify intensity range:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Specify intensity range:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    controls.lower = threshold_entry_attach(&controls, table,
                             row, args->lower, _("_Lower:"));
    g_signal_connect_swapped(controls.lower, "activate",
                             G_CALLBACK(threshold_lower_changed), &controls);
    row++;
    controls.upper = threshold_entry_attach(&controls, table, row,
                            args->upper, _("_Upper:"));
    g_signal_connect_swapped(controls.upper, "activate",
                            G_CALLBACK(threshold_upper_changed), &controls);
    row++;
    button = gtk_button_new_with_mnemonic(_("Set to _Full Range"));
    gtk_table_attach(table, button, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(threshold_set_to_full_range),
                             &controls);
    row++;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 20);
    label = gtk_label_new("Peak Positions:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Peak positions:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    GtkTreeViewColumn *column;
    GtkCellRenderer *renderer;
    GwyNullStore *store;
    store = gwy_null_store_new(2);
    tool->model = GTK_TREE_MODEL(store);
    tool->treeview = GTK_TREE_VIEW(gtk_tree_view_new_with_model(tool->model));
    gchar *XUnits, *YUnits, *ZUnits;
    XUnits = g_strdup_printf("<b>x</b> [%s]", controls.XY_Format->units);
    YUnits = g_strdup_printf("<b>y</b> [%s]", controls.XY_Format->units);
    ZUnits = g_strdup_printf("<b>value</b> [%s]", controls.Z_Format->units);
    guint i;
    for (i = 0; i < NCOLUMNS; i++) {
        column = gtk_tree_view_column_new();
        g_object_set_data(G_OBJECT(column), "id", GUINT_TO_POINTER(i));
        renderer = gtk_cell_renderer_text_new();
        g_object_set(renderer, "xalign", 1.0, NULL);
        gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(column), renderer, TRUE);
        gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(column), renderer,
                                            gwy_tool_level3_render_cell,
                                            &controls, NULL);
        label = gtk_label_new(NULL);
        switch (i)
        {
            case 0:
                gtk_label_set_markup(GTK_LABEL(label), "<b>n</b>");
                break;
            case 1:
                gtk_label_set_markup(GTK_LABEL(label), XUnits);
                break;
            case 2:
                gtk_label_set_markup(GTK_LABEL(label), YUnits);
                break;
            case 3:
                gtk_label_set_markup(GTK_LABEL(label), ZUnits);
                break;
        }
        gtk_tree_view_column_set_widget(column, label);
        gtk_widget_show(label);
        gtk_tree_view_append_column(tool->treeview, column);
    }
    gtk_table_attach(table, GTK_WIDGET(tool->treeview), 0, 3, row, row+1,
                        GTK_FILL, 0, 0, 0);
    row++;
    g_free(XUnits);
    g_free(YUnits);
    g_free(ZUnits);
    button = gtk_button_new_with_mnemonic(_("Clear Points"));
    gtk_table_attach(table, button, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(clear_points), &controls);
    row++;
    tool->radius = gtk_adjustment_new(tool->rpx, 0, 10, 1, 5, 0);
    gwy_table_attach_spinbutton(GTK_WIDGET(table), row, 
                        _("Peak search radius:"), "px", tool->radius);
    g_signal_connect_swapped(tool->radius, "value-changed",
                         G_CALLBACK(gwy_tool_level3_radius_changed), tool);
    row++;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 20);
    label = gtk_label_new("Specify HCP lattice constant:");
    gtk_label_set_markup(GTK_LABEL(label),
                        "<b>Specify HCP lattice constant:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    controls.lattice = threshold_entry_attach(&controls, table,
                        row, args->lattice, _("Lattice Constant:"));
    g_signal_connect_swapped(controls.lattice, "activate",
                         G_CALLBACK(threshold_lattice_changed), &controls);
    row++;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 20);
    label = gtk_label_new("Scale Factors:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Scale Factors:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    scale_entry_attach(&controls, table, row);
    row += 2;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 5);
    controls.warning = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(controls.warning), 0.5, 0.5);
    gtk_table_attach(table, controls.warning, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    preview(&controls);
    gtk_widget_show_all(dialog);
    do
    {
        response = gtk_dialog_run(GTK_DIALOG(dialog));
        switch (response)
        {
            case GTK_RESPONSE_CANCEL:
            case GTK_RESPONSE_DELETE_EVENT:
                gtk_widget_destroy(dialog);
            case GTK_RESPONSE_NONE:
                g_object_unref(controls.mydata);
                gwy_si_unit_value_format_free(controls.XY_Format);
                gwy_si_unit_value_format_free(controls.Z_Format);
                threshold_save_args(gwy_app_settings_get(), args, tool);
                return;
                break;
            case GTK_RESPONSE_OK:
                break;
            default:
                g_assert_not_reached();
                break;
        }
    } while (response != GTK_RESPONSE_OK);
    threshold_save_args(gwy_app_settings_get(), args, tool);
    if (gwy_selection_is_full(controls.selection))
        calibrate_do(&controls);
    else if (controls.args->Xscale > 0 && controls.args->Yscale > 0)
        calibrate_do(&controls);
    gtk_widget_destroy(dialog);
    g_object_unref(controls.mydata);
    gwy_si_unit_value_format_free(controls.original_XY_Format);
    gwy_si_unit_value_format_free(controls.XY_Format);
    gwy_si_unit_value_format_free(controls.Z_Format);
}

static void
threshold_set_to_range(ThresholdControls *controls,
                       gdouble lower, gdouble upper)
{
    threshold_format_value(controls, GTK_ENTRY(controls->lower), lower);
    gtk_widget_activate(controls->lower);
    threshold_format_value(controls, GTK_ENTRY(controls->upper), upper);
    gtk_widget_activate(controls->upper);
    preview(controls);
}

static void
threshold_set_to_full_range(ThresholdControls *controls)
{
    threshold_set_to_range(controls,
                   controls->ranges->min, controls->ranges->max);
}

static void
threshold_lower_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->lower));
    gdouble num =
        g_strtod(value, NULL) * controls->original_XY_Format->magnitude;
    if (num >= controls->ranges->min && num <= controls->ranges->max)
        controls->args->lower = num;
    else
    {
        if (num < controls->ranges->min)
            controls->args->lower = controls->ranges->min;
        else if (num > controls->ranges->max)
            controls->args->lower = controls->ranges->max; 
    }
    threshold_format_value(controls,
        GTK_ENTRY(controls->lower), controls->args->lower);
    preview(controls);
}

static void
threshold_upper_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->upper));
    gdouble num =
        g_strtod(value, NULL) * controls->original_XY_Format->magnitude;
    if (num >= controls->ranges->min && num <= controls->ranges->max)
        controls->args->upper = num;
    else
    {
        if (num < controls->ranges->min)
            controls->args->upper = controls->ranges->min;
        else if (num > controls->ranges->max)
            controls->args->upper = controls->ranges->max;
    }
    threshold_format_value(controls, GTK_ENTRY(controls->upper),
            controls->args->upper);
    preview(controls);
}

static void
threshold_lattice_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->lattice));
    gdouble num = g_strtod(value, NULL) * controls->original_XY_Format->magnitude;
    if (num > 0)
    {
        controls->args->lattice = num;
        calibrate_update_scales(controls);
    }
    else
        threshold_format_value(controls, GTK_ENTRY(controls->lattice),
                                                controls->args->lattice);
}

static void
preview(ThresholdControls *controls)
{
    GwyDataField *dfield;
    gdouble Xreal, Yreal;
    gdouble Xoff, Yoff;
    GwySIUnit *XY_Units;
    GwySIUnit *Z_Units;
    Xreal = gwy_data_field_get_xreal(controls->offt);
    Yreal = gwy_data_field_get_yreal(controls->offt);
    Xoff = gwy_data_field_get_xoffset(controls->offt);
    Yoff = gwy_data_field_get_yoffset(controls->offt);
    XY_Units = gwy_data_field_get_si_unit_xy(controls->offt);
    Z_Units = gwy_data_field_get_si_unit_z(controls->offt);
    ZoomMode zoom = controls->args->zoom_mode;
    dfield = GWY_DATA_FIELD(
            gwy_container_get_object_by_name(controls->mydata, "/0/data"));
    if (zoom != ZOOM_1)
    {
        gint Xres, Yres;
        Xres = gwy_data_field_get_xres(controls->offt);
        Yres = gwy_data_field_get_yres(controls->offt);
        guint width = (Xres/controls->args->zoom_mode) | 1;
        guint height = (Yres/controls->args->zoom_mode) | 1;
        GwyDataField *temp = gwy_data_field_area_extract(controls->offt,
                                          (Xres - width)/2, (Yres - height)/2,
                                          width, height);
        gwy_data_field_resample(temp, Xres, Yres, GWY_INTERPOLATION_BILINEAR);
        gwy_data_field_copy(temp, controls->disp_data, FALSE);
        g_object_unref(temp);
    }
    else
        gwy_data_field_copy(controls->offt, controls->disp_data, FALSE);
    gwy_data_field_set_xreal(controls->disp_data, Xreal/zoom);
    gwy_data_field_set_yreal(controls->disp_data, Yreal/zoom);
    gwy_data_field_set_xoffset(controls->disp_data, Xoff/zoom);
    gwy_data_field_set_yoffset(controls->disp_data, Yoff/zoom);
    gwy_data_field_set_si_unit_xy(controls->disp_data, XY_Units);
    gwy_data_field_set_si_unit_z(controls->disp_data, Z_Units);
    gwy_data_field_copy(controls->disp_data, dfield, FALSE);
    threshold_do(controls->args, dfield);
}

static void
reFind_Peaks(ThresholdControls *controls)
{
    int i, num = 0;
    for (i = 0; i < 4; i++)
    {
        double point[2];
        if (gwy_selection_get_object(controls->selection, i, point))
        {
            gdouble xoff, yoff;
            peak_find(controls, point, i);
            xoff = gwy_data_field_get_xoffset(controls->disp_data);
            yoff = gwy_data_field_get_yoffset(controls->disp_data);
            point[0] = controls->p[num][0] - xoff;
            point[1] = controls->p[num][1] - yoff;
            gwy_selection_set_object(controls->selection, i, point);
            num++;
        }
    }
    preview(controls);
}

static void
zoom_adjust_peaks(ThresholdControls *controls)
{
    int i, num = 0;
    for (i = 0; i < 4; i++)
    {
        double point[2];
        if (gwy_selection_get_object(controls->selection, i, point))
        {
            gdouble xoff, yoff, multiplier;
            multiplier = 1.0/controls->args->zoom_mode;
            xoff = gwy_data_field_get_xoffset(controls->dfield) * multiplier;
            yoff = gwy_data_field_get_yoffset(controls->dfield) * multiplier;
            point[0] = controls->p[num][0] - xoff;
            point[1] = controls->p[num][1] - yoff;
            gwy_selection_set_object(controls->selection, i, point);
            num++;
        }
    }
    reFind_Peaks(controls);
}

static void
threshold_do(const ThresholdArgs *args, GwyDataField *dfield)
{
    gdouble lower = MIN(args->lower, args->upper);
    gdouble upper = MAX(args->lower, args->upper);
    gwy_data_field_clamp(dfield, lower, upper);
    gwy_data_field_data_changed(dfield);
}

static void
clear_points(ThresholdControls *controls)
{
    gwy_selection_clear(controls->selection);
    calibrate_update_scales(controls);
}

static const gchar lower_key[] = "/module/calibrate_hcp/lower";
static const gchar upper_key[] = "/module/calibrate_hcp/upper";
static const gchar lattice_key[] = "/module/calibrate_hcp/lattice";
static const gchar radius_key[] = "/module/calibrate_hcp/radius";

static void
threshold_load_args(GwyContainer *settings, 
                    ThresholdArgs *args, GwyToolLevel3 *tool)
{
    *args = threshold_defaults;
    gwy_container_gis_double_by_name(settings, lower_key, &args->lower);
    gwy_container_gis_double_by_name(settings, upper_key, &args->upper);
    gwy_container_gis_double_by_name(settings, lattice_key, &args->lattice);
    gwy_container_gis_int32_by_name(settings, radius_key, &(tool->rpx));
}

static void
threshold_save_args(GwyContainer *settings,
                    ThresholdArgs *args, GwyToolLevel3 *tool)
{
    gwy_container_set_double_by_name(settings, lower_key, args->lower);
    gwy_container_set_double_by_name(settings, upper_key, args->upper);
    gwy_container_set_double_by_name(settings, lattice_key, args->lattice);
    gwy_container_set_int32_by_name(settings, radius_key, tool->rpx);
}

static void
gwy_tool_level3_render_cell(GtkCellLayout *layout,
            GtkCellRenderer *renderer, GtkTreeModel *model,
            GtkTreeIter *iter, gpointer user_data)
{
    ThresholdControls *controls = (ThresholdControls*)user_data;
    const GwySIValueFormat *vf;
    gchar buf[32];
    gdouble point[2];
    gdouble val;
    guint idx, id;
    id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(layout), "id"));
    gtk_tree_model_get(model, iter, 0, &idx, -1);
    if (id == COLUMN_I)
    {
        g_snprintf(buf, sizeof(buf), "%d", idx + 1);
        g_object_set(renderer, "text", buf, NULL);
        return;
    }
    if (!controls->selection ||
            !gwy_selection_get_object(controls->selection, idx, point))
    {
        g_object_set(renderer, "text", "", NULL);
        return;
    }
    switch (id)
    {
        case COLUMN_X:
            peak_find(controls, point, idx);
            vf = controls->XY_Format;
            val = controls->p[idx][0];
            break;
        case COLUMN_Y:
            vf = controls->XY_Format;
            val = controls->p[idx][1];
            break;
        case COLUMN_Z:
            vf = controls->Z_Format;
            val = controls->p[idx][2];
            break;
        default:
            g_return_if_reached();
            break;
    }
    if (vf)
        g_snprintf(buf, sizeof(buf), "%.*f",
            vf->precision, val/vf->magnitude);
    else
        g_snprintf(buf, sizeof(buf), "%.3g", val);
    g_object_set(renderer, "text", buf, NULL);
    calibrate_update_scales(controls);
}

static void
calibrate_update_scales(ThresholdControls *controls)
{
    if (gwy_selection_is_full(controls->selection))
    {
        calibration_get_factors(controls);
        gtk_entry_set_text(GTK_ENTRY(controls->xscale),
            g_strdup_printf("%f", controls->args->Xscale));
        gtk_entry_set_text(GTK_ENTRY(controls->yscale),
            g_strdup_printf("%f", controls->args->Yscale));
    }
    else
    {
        gchar *s1 = g_strdup_printf("%f", controls->args->Xscale);
        gchar *s2 = g_strdup_printf("%f", controls->args->Yscale);
        controls->args->Xscale = 0.0;
        controls->args->Yscale = 0.0;
        gtk_entry_set_text(GTK_ENTRY(controls->xscale), s1);
        gtk_entry_set_text(GTK_ENTRY(controls->yscale), s2);
        controls->args->Xwarning = FALSE;
        controls->args->Ywarning = FALSE;
        check_warnings(controls);
        g_free(s1);
        g_free(s2);
    }
}

static void
peak_find(ThresholdControls *controls, gdouble *point, guint idx)
{
    GwyDataField *dfield = controls->disp_data;
    gint i, j, low_i, high_i, low_j, high_j;
    gdouble temp_i, temp_j, temp_z;
    gint col = gwy_data_field_rtoj(dfield, point[0]);
    gint row = gwy_data_field_rtoi(dfield, point[1]);
    temp_i = col;
    temp_j = row;    
    temp_z = gwy_data_field_get_val(dfield, col, row);
    gint32 r = controls->tool->rpx;
    gint Xres, Yres;
    Xres = gwy_data_field_get_xres(dfield);
    Yres = gwy_data_field_get_yres(dfield);
    low_i = col - r;
    high_i = col + r;
    if (low_i < 0)
        low_i = 0;
    if (high_i > Xres)
        high_i = Xres;
    low_j = row - r;
    high_j = row + r;
    if (low_j < 0)
        low_j = 0;
    if (high_j > Yres)
        high_j = Yres;
    for (i = low_i; i < high_i; i++)
    {
        for (j = low_j; j < high_j; j++)
        {
            if (gwy_data_field_get_val(dfield, i, j) > temp_z)
            {
                temp_i = i;
                temp_j = j;
                temp_z = gwy_data_field_get_val(dfield, i, j);
            }
        }
    }
    controls->p[idx][0] = gwy_data_field_jtor(dfield, temp_i)
                + gwy_data_field_get_xoffset(dfield);
    controls->p[idx][1] = gwy_data_field_itor(dfield, temp_j)
                + gwy_data_field_get_yoffset(dfield);
    controls->p[idx][2] = temp_z;
    if ((row - temp_j) != 0 || (col - temp_i) != 0)
    {
        point[0] = gwy_data_field_jtor(dfield, temp_i);
        point[1] = gwy_data_field_itor(dfield, temp_j);
        gwy_selection_set_object(controls->selection, idx, point);
    }
}

static void
gwy_tool_level3_radius_changed(GwyToolLevel3 *tool)
{
    tool->rpx = gwy_adjustment_get_int(tool->radius);
    guint i;
    GwyNullStore *store = GWY_NULL_STORE(tool->model);
    for (i = 0; i < 2; i++)
        gwy_null_store_row_changed(store, i);
}

static void
perform_fft(GwyDataField *dfield, GwyContainer *data)
{    
    GwyDataField *raout, *ipout;
    raout = gwy_data_field_new_alike(dfield, FALSE);
    ipout = gwy_data_field_new_alike(dfield, FALSE);
    gwy_data_field_2dfft(dfield, NULL, raout, ipout,
                         GWY_WINDOWING_HANN,
                         GWY_TRANSFORM_DIRECTION_FORWARD,
                         GWY_INTERPOLATION_LINEAR, FALSE, 1);
    set_dfield_modulus(raout, ipout, dfield);
    fft_postprocess(dfield);
    gchar *key;
    key = g_strdup_printf("/%i/base/palette", 0);
    gwy_container_set_string_by_name(data, key, g_strdup("Gray"));
    g_free(key);
    key = g_strdup_printf("/%i/base/range-type", 0);
    gwy_container_set_enum_by_name(data, key, GWY_LAYER_BASIC_RANGE_ADAPT);
    g_free(key);
    g_object_unref(raout);
    g_object_unref(ipout);
}

static void
set_dfield_modulus(GwyDataField *re, GwyDataField *im, GwyDataField *target)
{
    const gdouble *datare, *dataim;
    gdouble *data;
    gint xres, yres, i;
    xres = gwy_data_field_get_xres(re);
    yres = gwy_data_field_get_yres(re);
    datare = gwy_data_field_get_data_const(re);
    dataim = gwy_data_field_get_data_const(im);
    data = gwy_data_field_get_data(target);
    for (i = xres*yres; i; i--, datare++, dataim++, data++)
        *data = hypot(*datare, *dataim);
}

static void
fft_postprocess(GwyDataField *dfield)
{
    gint res;
    gdouble r;
    gwy_data_field_2dfft_humanize(dfield);
    GwySIUnit *xyunit;
    xyunit = gwy_data_field_get_si_unit_xy(dfield);
    gwy_si_unit_power(xyunit, -1, xyunit);
    gwy_data_field_set_xreal(dfield, 1.0/gwy_data_field_get_xmeasure(dfield));
    gwy_data_field_set_yreal(dfield, 1.0/gwy_data_field_get_ymeasure(dfield));
    res = gwy_data_field_get_xres(dfield);
    r = res / 2.0;
    gwy_data_field_set_xoffset(dfield, -gwy_data_field_jtor(dfield, r));
    res = gwy_data_field_get_yres(dfield);
    r = res / 2.0;
    gwy_data_field_set_yoffset(dfield, -gwy_data_field_itor(dfield, r));
    gdouble dmin, dmax;
    gwy_data_field_get_min_max(dfield, &dmin, &dmax);
    gwy_data_field_add(dfield, -dmin);
}

static void
selection_changed(ThresholdControls *controls)
{
    guint i;
    GwyNullStore *store = GWY_NULL_STORE(controls->tool->model);
    for (i = 0; i < 2; i++)
        gwy_null_store_row_changed(store, i);
}

static void
calibration_get_factors(ThresholdControls *controls)
{

    gdouble x1, x2, y1, y2, R, xcorr, ycorr;
    gdouble x1_2, x2_2, y1_2, y2_2;
    x1 = controls->p[0][0];
    y1 = controls->p[0][1];
    x2 = controls->p[1][0];
    y2 = controls->p[1][1];
    R = 2 / (sqrt(3) * controls->args->lattice);
    x1_2 = x1 * x1;
    y1_2 = y1 * y1;
    x2_2 = x2 * x2;
    y2_2 = y2 * y2;
    ycorr = R * sqrt((x1_2 - x2_2) / (x1_2 * y2_2 - x2_2 * y1_2));
    xcorr = sqrt((R * R - (ycorr * ycorr * y1_2)) / x1_2);
    controls->args->Xscale = 1 / xcorr;
    controls->args->Yscale = 1 / ycorr;
    controls->args->Xwarning = FALSE;
    controls->args->Ywarning = FALSE;
    if (x1_2 == x2_2)
        controls->args->Xwarning = TRUE;
    if (y1_2 == y2_2)
        controls->args->Ywarning = TRUE;
    if (x1_2 == 0)
        controls->args->Xwarning = TRUE;
    if (controls->args->Xscale != controls->args->Xscale)
        controls->args->Xwarning = TRUE;
    if (controls->args->Yscale != controls->args->Yscale)
        controls->args->Ywarning = TRUE;
    if (x1_2 * y2_2 == x2_2 * y1_2)
    {
        controls->args->Xwarning = TRUE;
        controls->args->Ywarning = TRUE;
    }
    check_warnings(controls);
}

static void
check_warnings(ThresholdControls *controls)
{
    if (controls->args->Xwarning)
        gtk_label_set_markup(GTK_LABEL(controls->xwarning),
            "<span foreground=\"red\"><b>X</b></span>");
    else
        gtk_label_set_markup(GTK_LABEL(controls->xwarning), "");
    if (controls->args->Ywarning)
        gtk_label_set_markup(GTK_LABEL(controls->ywarning),
            "<span foreground=\"red\"><b>X</b></span>");
    else
        gtk_label_set_markup(GTK_LABEL(controls->ywarning), "");
    if (controls->args->Xwarning || controls->args->Ywarning)
        gtk_label_set_markup(GTK_LABEL(controls->warning),
            "<span foreground=\"red\"><b>Warning!</b></span>");
    else
        gtk_label_set_markup(GTK_LABEL(controls->warning), "");
}

static void
calibrate_do(ThresholdControls *controls)
{
    gint oldXres = gwy_data_field_get_xres(controls->ofield);
    gint oldYres = gwy_data_field_get_yres(controls->ofield);
    gint newXres = GWY_ROUND(oldXres);
    gint newYres = GWY_ROUND(oldYres *
        controls->args->Yscale / controls->args->Xscale);
    GwyDataField *newDataField = gwy_data_field_new_resampled
        (controls->ofield, newXres, newYres, GWY_INTERPOLATION_LINEAR);
    gdouble oldXreal = gwy_data_field_get_xreal(controls->ofield);
    gdouble oldYreal = gwy_data_field_get_yreal(controls->ofield);
    gdouble newXreal = oldXreal * controls->args->Xscale;
    gdouble newYreal = oldYreal * controls->args->Yscale;
    gwy_data_field_set_xreal(newDataField, newXreal);
    gwy_data_field_set_yreal(newDataField, newYreal);
    calibrate_create_output(controls->container, newDataField, controls);
}

static void
calibrate_create_output(GwyContainer *data,
    GwyDataField *dfield, ThresholdControls *controls)
{
    gint id, newid;
    const guchar *title;
    GwyContainer *meta;
    gwy_app_data_browser_get_current(GWY_APP_DATA_FIELD_ID, &id, 0);
    GQuark Qmeta = g_quark_from_string(g_strdup_printf("/%i/meta", id));
    if (gwy_container_contains(data, Qmeta))
        meta = gwy_container_duplicate(gwy_container_get_object(data, Qmeta));
    else
        meta = gwy_container_new();
    title = gwy_container_get_string(data,
            g_quark_try_string(g_strdup_printf("/%i/data/title", id)));
    gwy_container_set_string_by_name(meta, "Source Title", title);
    gwy_container_set_string_by_name(meta, "X Scaling Factor",
            (const guchar *)g_strdup_printf("%.5f", controls->args->Xscale));
    gwy_container_set_string_by_name(meta, "Y Scaling Factor",
            (const guchar *)g_strdup_printf("%.5f", controls->args->Yscale));
    newid = gwy_app_data_browser_add_data_field(dfield, data, TRUE);
    gwy_container_set_object_by_name(data,
            g_strdup_printf("/%i/meta", newid), meta);
    gwy_app_set_data_field_title(data, newid, _("Calibrated"));
    gwy_app_channel_log_add(data, controls->id,
            newid, "proc::calibrate_hcp", NULL);
    g_object_unref(dfield);
}

static void
zoom_mode_changed(GtkToggleButton *button, ThresholdControls *controls)
{
    controls->args->zoom_mode =
        gwy_radio_buttons_get_current(controls->zoom_mode_radios);
    preview(controls);
    zoom_adjust_peaks(controls);
}

static void
radio_buttons_attach_to_table(GSList *group, 
                GtkTable *table, gint row)
{
    g_return_val_if_fail(GTK_IS_TABLE(table), row);
    while (group)
    {
        gtk_table_attach(table, GTK_WIDGET(group->data),
                         0, 1, row, row + 1,
                         GTK_EXPAND | GTK_FILL, 0, 0, 0);
        group = g_slist_next(group);
        gtk_table_attach(table, GTK_WIDGET(group->data),
                         2, 3, row, row + 1,
                         GTK_EXPAND | GTK_FILL, 0, 0, 0);
        row++;
        group = g_slist_next(group);
    }
}

static void
xscale_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->xscale));
    gdouble num = g_strtod(value, NULL);
    if (num != controls->args->Xscale && num > 0)
    {
        controls->args->Xscale = num;
    }
    gchar *s = g_strdup_printf("%0.1f", controls->args->Xscale);
    gtk_entry_set_text(GTK_ENTRY(controls->xscale), s);
    g_free(s);
}

static void
yscale_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->yscale));
    gdouble num = g_strtod(value, NULL);
    if (num != controls->args->Yscale && num > 0)
    {
        controls->args->Yscale = num;
    }
    gchar *s = g_strdup_printf("%0.1f", controls->args->Yscale);
    gtk_entry_set_text(GTK_ENTRY(controls->yscale), s);
    g_free(s);
}
