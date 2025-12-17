"""
Helper functions for managing tooltip configurations in plots.
"""

from multiqc import config


def should_disable_tooltips():
    """
    Check if tooltips should be disabled based on the global configuration.
    
    Returns:
        bool: True if tooltips should be disabled, False otherwise
    """
    return getattr(config, 'disable_plot_tooltips', False)


def apply_tooltip_config(plot_config):
    """
    Apply tooltip configuration to a plot config dictionary.
    
    For plots that use tt_label or tt_decimals, this will remove or modify
    tooltip-related configurations based on the global disable_plot_tooltips setting.
    
    Args:
        plot_config (dict): The plot configuration dictionary to modify
        
    Returns:
        dict: The modified plot configuration dictionary
    """
    if should_disable_tooltips():
        # Remove tooltip-related keys if they exist
        plot_config.pop('tt_decimals', None)
        plot_config.pop('tt_suffix', None)
        plot_config.pop('tt_percentages', None)
        
        # Set tt_label to False to disable tooltips
        # This works for most plot types in MultiQC/pmultiqc
        plot_config['tt_label'] = False
    
    return plot_config
