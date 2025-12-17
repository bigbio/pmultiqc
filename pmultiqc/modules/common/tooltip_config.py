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
    # Use getattr with fallback for safety in case config attribute doesn't exist yet
    return getattr(config, 'disable_plot_tooltips', False)


def apply_tooltip_config(plot_config):
    """
    Apply tooltip configuration to a plot config dictionary or data label dictionary.
    
    For plots that use tt_label or tt_decimals, this will remove or modify
    tooltip-related configurations based on the global disable_plot_tooltips setting.
    
    This function can be applied to:
    - Plot configuration dictionaries (pconfig)
    - Data label dictionaries within plot configs
    
    Args:
        plot_config (dict): The plot configuration dictionary or data label to modify
        
    Returns:
        dict: The modified dictionary
    """
    if should_disable_tooltips():
        # Remove tooltip-related keys if they exist
        plot_config.pop('tt_decimals', None)
        plot_config.pop('tt_suffix', None)
        plot_config.pop('tt_percentages', None)
        
        # Set tt_label to False to disable tooltips
        # This works for both plot configs and data labels in MultiQC/pmultiqc
        plot_config['tt_label'] = False
    
    return plot_config
